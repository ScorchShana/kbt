#ifndef FASTA_READER_HEADER
#define FASTA_READER_HEADER

#include "definition.h"
#include "RingMemoryPool.h"
#include "SPSCRingQueue.h"
#include "SpinBackoff.h"

#include <cassert>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <utility>
#include <iostream>
#include <cstdlib>
#include <zlib.h>
#include <thread>
#include <atomic>
#include <algorithm>

class FASTAReaderThreadPool
{
    static constexpr uint64_t kNoNewlineInBlock = static_cast<uint64_t>(-1);
    static constexpr uint64_t GZ_CHUNK_SIZE = 1 * 1024 * 1024;

    int k_;
    uint64_t base_chunk_size_;
    RingMemoryPool<READER_PARSER_RING_MEMORY_POOL_CAPACITY>* ring_pool_ptr_;

    std::vector<std::string> gz_files_;
    std::vector<std::string> plain_files_;
    uint32_t reader_count_;
    std::vector<std::unique_ptr<std::thread>> threads_;
    std::vector<std::string> assignments;

    alignas(CACHE_LINE_SIZE) std::atomic<uint64_t> file_index_{ 0 };

    inline static thread_local SpinBackoff<> enqueue_backoff;
    inline static thread_local SpinBackoff<> dequeue_backoff;

public:

#ifdef TEST_MODE
    std::atomic<uint64_t> total_dequeue_spin_time = 0;
    std::atomic<uint64_t> total_enqueue_spin_time = 0;
    inline thread_local static uint64_t dequeue_spin_time = 0;
    inline thread_local static uint64_t enqueue_spin_time = 0;
#endif

    explicit FASTAReaderThreadPool(
        const std::vector<std::string>& filenames,
        uint64_t chunk_size,
        uint32_t reader_count,
        RingMemoryPool<READER_PARSER_RING_MEMORY_POOL_CAPACITY>* pool_ptr)
        : k_(k_length), base_chunk_size_(chunk_size), ring_pool_ptr_(pool_ptr), reader_count_(reader_count)
    {
        assert(ring_pool_ptr_ != nullptr);
        assert(ring_pool_ptr_->blockSize() > 1024);
        if (filenames.empty()) { std::cerr << "No input files" << std::endl; throw std::runtime_error("input I/O failed"); }

        classify_files(filenames);
        sort_files_by_size();

        threads_.reserve(reader_count_);
    }

    ~FASTAReaderThreadPool()
    {
        RunState::fail_on_unwind();
        join();
#ifdef TEST_MODE
        std::cout << "FASTAReaderThreadPool: total dequeue spin time: " << total_dequeue_spin_time.load(std::memory_order_relaxed) << std::endl;
        std::cout << "FASTAReaderThreadPool: total enqueue spin time: " << total_enqueue_spin_time.load(std::memory_order_relaxed) << std::endl;
#endif
    }

    void start()
    {
        file_index_.store(0, std::memory_order_relaxed);
        assignments.clear();
        for (size_t i = 0; i < gz_files_.size(); ++i)
            assignments.push_back(gz_files_[i]);
        for (size_t i = 0; i < plain_files_.size(); ++i)
            assignments.push_back(plain_files_[i]);

        for (uint32_t i = 0; i < reader_count_; ++i)
        {
            threads_.push_back(std::make_unique<std::thread>([this]() {
                reader_worker();
                }));
        }
    }

    void join()
    {
        for (auto& t : threads_)
        {
            if (t->joinable()) t->join();
        }
    }

private:
    void classify_files(const std::vector<std::string>& filenames)
    {
        for (const auto& f : filenames)
        {
            if (f.size() >= 3 && f.compare(f.size() - 3, 3, ".gz") == 0)
                gz_files_.push_back(f);
            else
                plain_files_.push_back(f);
        }
    }

    void sort_files_by_size()
    {
        struct file_info {
            std::string filename;
            uint64_t size;
        };
        std::vector<file_info> infos;

        for (const auto& gz_filename : gz_files_)
        {
            struct stat st;
            uint64_t file_size = 0;
            if (stat(gz_filename.c_str(), &st) == -1) file_size = 0;
            else file_size = static_cast<uint64_t>(st.st_size);
            infos.push_back({ gz_filename, file_size });
        }

        std::sort(infos.begin(), infos.end(), [](const file_info& a, const file_info& b) {
            return a.size > b.size;
            });
        gz_files_.clear();
        for (const auto& info : infos)
        {
            gz_files_.push_back(info.filename);
        }

        infos.clear();
        for (const auto& plain_filename : plain_files_)
        {
            struct stat st;
            uint64_t file_size = 0;
            if (stat(plain_filename.c_str(), &st) == -1) file_size = 0;
            else file_size = static_cast<uint64_t>(st.st_size);
            infos.push_back({ plain_filename, file_size });
        }

        std::sort(infos.begin(), infos.end(), [](const file_info& a, const file_info& b) {
            return a.size > b.size;
            });
        plain_files_.clear();
        for (const auto& info : infos)
        {
            plain_files_.push_back(info.filename);
        }


    }

    inline void acquire_block(char*& block_ptr, uint64_t& write_size, uint64_t& last_newline_pos,
        bool& has_block, size_t& left_buffer_size_, char* left_buffer_)
    {
        if (ring_pool_ptr_->producer_try_dequeue(block_ptr))
        {
            dequeue_backoff.reset();
        }
        else
        {
            while (!ring_pool_ptr_->producer_try_dequeue(block_ptr))
            {
                dequeue_backoff.backoff();
#ifdef TEST_MODE
                ++dequeue_spin_time;
#endif
            }
            dequeue_backoff.decay();
        }


        has_block = true; write_size = 0; last_newline_pos = kNoNewlineInBlock;
        if (left_buffer_size_ > 0)
        {
            std::memcpy(block_ptr, left_buffer_, left_buffer_size_);
            write_size = left_buffer_size_;
            left_buffer_size_ = 0;
        }
    }

    inline void publish_current_block(char*& block_ptr, uint64_t& write_size,
        uint64_t& last_newline_pos, bool& has_block)
    {
        if (!has_block) return;

        if (write_size > 0) {
            if (ring_pool_ptr_->producer_try_enqueue({ block_ptr, write_size }))
            {
                enqueue_backoff.reset();
            }
            else
            {
                while (!ring_pool_ptr_->producer_try_enqueue({ block_ptr, write_size }))
                {
                    enqueue_backoff.backoff();
#ifdef TEST_MODE
                    ++enqueue_spin_time;
#endif
                }
                enqueue_backoff.decay();
            }
        }
        else ring_pool_ptr_->consumer_enqueue(block_ptr);

        has_block = false;
        block_ptr = nullptr;
        write_size = 0;
        last_newline_pos = kNoNewlineInBlock;
    }

    inline void store_overlap_from_block_end(const char* block_ptr, uint64_t write_size,
        uint64_t last_newline_pos, uint64_t overlap,
        size_t& left_buffer_size_, char* left_buffer_)
    {
        left_buffer_size_ = 0;
        if (overlap == 0 || write_size == 0) return;
        const uint64_t keep = (write_size < overlap) ? write_size : overlap;
        const uint64_t start = write_size - keep;
        if (last_newline_pos != kNoNewlineInBlock && last_newline_pos >= start) {
            const uint64_t kp = write_size - (last_newline_pos + 1);
            if (kp > 0) { std::memcpy(left_buffer_, block_ptr + last_newline_pos + 1, kp); left_buffer_size_ = static_cast<size_t>(kp); }
            return;
        }
        std::memcpy(left_buffer_, block_ptr + start, keep);
        left_buffer_size_ = static_cast<size_t>(keep);
    }

    void reader_worker()
    {
        run_worker("reader", [&] {

        assert(ring_pool_ptr_ != nullptr);
        const uint64_t block_size = ring_pool_ptr_->blockSize();

        char left_buffer_[128];
        size_t left_buffer_size_ = 0;
        const uint64_t overlap = (k_ > 1) ? static_cast<uint64_t>(k_ - 1) : 0;
        const uint64_t worker_buffer_size = std::max(base_chunk_size_, GZ_CHUNK_SIZE);
        std::vector<char> read_buf(worker_buffer_size);

        uint64_t cur_file_index = file_index_.fetch_add(1, std::memory_order_relaxed);
        for (; cur_file_index < assignments.size(); cur_file_index = file_index_.fetch_add(1, std::memory_order_relaxed))
        {
            const std::string& file = assignments[cur_file_index];
            const bool is_gz = (file.size() >= 3 && file.compare(file.size() - 3, 3, ".gz") == 0);
            const uint64_t effective_chunk_size = is_gz ? GZ_CHUNK_SIZE : base_chunk_size_;

            int fd = -1;
            gzFile gzfile = nullptr;

            if (is_gz)
            {
                gzfile = gzopen(file.c_str(), "rb");
                if (gzfile == nullptr) { std::cerr << "Failed to open gzip: " << file << std::endl; throw std::runtime_error("input I/O failed"); }
                gzbuffer(gzfile, GZ_CHUNK_SIZE / 4);
            }
            else
            {
                fd = ::open(file.c_str(), O_RDONLY);
                if (fd == -1) { std::cerr << "Failed to open: " << file << std::endl; throw std::runtime_error("input I/O failed"); }
                posix_fadvise(fd, 0, 0, POSIX_FADV_SEQUENTIAL);
            }

            bool line_start = true, in_header = false, have_record = false;
            left_buffer_size_ = 0;
            char* block_ptr = nullptr;
            uint64_t write_size = 0;
            bool has_block = false;
            uint64_t last_newline_pos = kNoNewlineInBlock;

            uint64_t input_pos = 0, input_size = 0;
            bool eof = false;
            char* input_begin = read_buf.data();

            try {
                while (true)
                {
                    if (input_pos >= input_size && !eof)
                    {
                        ssize_t bytes_read;
                        if (is_gz)
                        {
                            bytes_read = gzread(gzfile, read_buf.data(), static_cast<unsigned int>(effective_chunk_size));
                            input_begin = read_buf.data();
                        }
                        else
                        {
                            bytes_read = ::read(fd, read_buf.data(), effective_chunk_size);
                            input_begin = read_buf.data();
                        }


                        if (bytes_read < 0) [[unlikely]]
                        {
                            std::cerr << "Failed to read fasta data: " << file << std::endl;
                            throw std::runtime_error("input I/O failed");
                        }
                        if (bytes_read <= 0)
                        {
                            eof = true;
                        }
                        else
                        {
                            input_pos = 0;
                            input_size = static_cast<uint64_t>(bytes_read);
                        }
                    }

                    if (!has_block)
                        acquire_block(block_ptr, write_size, last_newline_pos, has_block, left_buffer_size_, left_buffer_);

                    if (input_pos >= input_size)
                    {
                        if (eof) {
                            if (have_record) {
                                if (write_size == block_size) {
                                    left_buffer_size_ = 0;
                                    publish_current_block(block_ptr, write_size, last_newline_pos, has_block);
                                    acquire_block(block_ptr, write_size, last_newline_pos, has_block, left_buffer_size_, left_buffer_);
                                }
                                block_ptr[write_size++] = '\n';
                            }
                            publish_current_block(block_ptr, write_size, last_newline_pos, has_block);
                            break;
                        }
                        continue;
                    }

                    if (line_start && input_begin[input_pos] == '>')
                    {
                        if (have_record)
                        {
                            // A record boundary needs no overlap, even at a full block.
                            if (write_size == block_size) {
                                left_buffer_size_ = 0;
                                publish_current_block(block_ptr, write_size, last_newline_pos, has_block);
                                acquire_block(block_ptr, write_size, last_newline_pos, has_block, left_buffer_size_, left_buffer_);
                            }
                            last_newline_pos = write_size;
                            block_ptr[write_size++] = '\n';
                        }
                        have_record = true;
                        in_header = true;
                    }
                    const char* cur = input_begin + input_pos;
                    const char* nl = static_cast<const char*>(std::memchr(cur, '\n', input_size - input_pos));
                    const uint64_t line_bytes = nl ? static_cast<uint64_t>(nl - cur) : input_size - input_pos;
                    if (!in_header)
                    {
                        // CR can be the final byte of an input chunk; do not emit it.
                        const uint64_t seq_len = line_bytes - (line_bytes > 0 && cur[line_bytes - 1] == '\r');
                        uint64_t copied = 0;
                        while (copied < seq_len)
                        {
                            if (write_size == block_size) {
                                store_overlap_from_block_end(block_ptr, write_size, last_newline_pos, overlap, left_buffer_size_, left_buffer_);
                                publish_current_block(block_ptr, write_size, last_newline_pos, has_block);
                                acquire_block(block_ptr, write_size, last_newline_pos, has_block, left_buffer_size_, left_buffer_);
                            }
                            const uint64_t take = std::min(seq_len - copied, block_size - write_size);
                            std::memcpy(block_ptr + write_size, cur + copied, take);
                            copied += take;
                            write_size += take;
                        }
                    }
                    input_pos += line_bytes + (nl != nullptr);
                    line_start = nl != nullptr;
                    if (line_start) in_header = false;
                }

            } catch (...) {
                if (gzfile) gzclose(gzfile);
                if (fd != -1) ::close(fd);
                throw;
            }

            if (is_gz) gzclose(gzfile);
            else ::close(fd);

            std::cout << "FASTAReader: completed " << file << std::endl;
        }

        ring_pool_ptr_->producer_set_finished();

#ifdef TEST_MODE
        total_dequeue_spin_time.fetch_add(dequeue_spin_time, std::memory_order_relaxed);
        total_enqueue_spin_time.fetch_add(enqueue_spin_time, std::memory_order_relaxed);
#endif

        });
    }
};

#endif
