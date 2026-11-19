#include "../src/FASTAPreReader.h"
#include <filesystem>
#include <fstream>
#include <future>
#include <map>
#include <random>
#include <string_view>

namespace {
void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}
using Counts = std::map<std::string, uint64_t>;
void count(std::string_view text, Counts& counts) {
    size_t valid = 0;
    for (size_t i = 0; i < text.size(); ++i) {
        char c = text[i];
        if (c == 'A' || c == 'C' || c == 'G' || c == 'T') ++valid;
        else valid = 0;
        if (valid >= k_length) ++counts[std::string(text.substr(i + 1 - k_length, k_length))];
    }
}
void check(const std::filesystem::path& dir, uint32_t, uint64_t chunk, uint64_t block) {
    std::mt19937 rng(73);
    std::vector<std::string> sequences{"", "A", "ACGT", std::string(block, 'A'), std::string(block * 3 + 17, 'C')};
    std::string random;
    for (int i = 0; i < 14000; ++i) random += "ACGTN"[rng() % 5];
    sequences.push_back(random);
    Counts expected, actual;
    std::vector<std::string> files;
    uint64_t plain_bytes = 0, gz_bytes = 0;
    for (int file = 0; file < 4; ++file) {
        std::string input;
        if (file < 2) {
            for (size_t i = 0; i < sequences.size(); ++i) {
                input += ">" + std::string(i == 0 ? (file == 1 ? 1024 * 1024 - 2 : 4094) : 7, 'h') + "\r\n";
                const auto& seq = sequences[i];
                const size_t width = file == 0 ? 31 : seq.size() + 1;
                for (size_t pos = 0; pos < seq.size(); pos += width)
                    input += seq.substr(pos, width) + "\r\n\n";
                count(seq, expected);
            }
            input += ">tail\nACGTACGTACGT"; // no final newline
            count("ACGTACGTACGT", expected);
        }
        bool gz = file % 2;
        auto path = dir / (std::to_string(file) + (gz ? ".fa.gz" : ".fa"));
        files.push_back(path.string());
        if (gz) {
            auto out = gzopen(path.c_str(), "wb");
            require(out != nullptr, "gzip fixture open");
            require(gzwrite(out, input.data(), input.size()) == static_cast<int>(input.size()), "gzip fixture write");
            require(gzclose(out) == Z_OK, "gzip fixture close");
            gz_bytes += input.size();
        } else {
            std::ofstream out(path, std::ios::binary);
            out << input;
            require(out.good(), "plain fixture write");
            plain_bytes += input.size();
        }
    }
    RunState state;
    ActiveRunScope active(state);
    RingMemoryPool<READER_PARSER_RING_MEMORY_POOL_CAPACITY> pool(block, 1);
    FASTAPreReader reader(files, chunk, &pool);
    auto worker = std::async(std::launch::async, [&] { run_worker("test pre-reader", [&] { reader.pre_read(); }); });
    content_type content;
    while (!pool.producer_finished()) {
        state.check();
        if (pool.consumer_try_dequeue(content)) {
            count({content.data, content.length}, actual);
            pool.consumer_enqueue(content.data);
        } else std::this_thread::yield();
    }
    worker.get();
    state.check();
    while (pool.consumer_try_dequeue(content)) {
        count({content.data, content.length}, actual);
        pool.consumer_enqueue(content.data);
    }
    require(actual == expected, "sequence boundary/overlap counts differ");
    require(reader.get_sampled_plain_bytes() == plain_bytes, "plain sampled bytes");
    require(reader.get_sampled_gz_decompressed_bytes() == gz_bytes, "gzip sampled bytes");
    require(reader.get_sampled_gz_compressed_bytes() > 0, "compressed sampled bytes");
}
void check_sampling_limit(const std::filesystem::path& dir) {
    const auto path = dir / "large.fa";
    constexpr uint64_t mib = 1024 * 1024;
    const std::string chunk(mib, 'A');
    {
        std::ofstream out(path, std::ios::binary);
        out << ">x\n";
        for (int i = 0; i < 257; ++i) out << chunk;
        require(out.good(), "sampling fixture write");
    }
    for (uint32_t files : {1U, 2U}) {
        RunState state;
        ActiveRunScope active(state);
        k_length = 31;
        RingMemoryPool<READER_PARSER_RING_MEMORY_POOL_CAPACITY> pool(4096, 1);
        FASTAPreReader reader(std::vector<std::string>(files, path.string()), mib, &pool);
        auto worker = std::async(std::launch::async, [&] { run_worker("sample limit", [&] { reader.pre_read(); }); });
        uint64_t kmers = 0;
        content_type content;
        while (!pool.producer_finished()) {
            state.check();
            if (!pool.consumer_try_dequeue(content)) { std::this_thread::yield(); continue; }
            uint64_t bases = content.length;
            if (bases && content.data[bases - 1] == '\n') --bases;
            if (bases >= k_length) kmers += bases - k_length + 1;
            pool.consumer_enqueue(content.data);
        }
        worker.get();
        state.check();
        while (pool.consumer_try_dequeue(content)) {
            uint64_t bases = content.length;
            if (bases && content.data[bases - 1] == '\n') --bases;
            if (bases >= k_length) kmers += bases - k_length + 1;
            pool.consumer_enqueue(content.data);
        }
        const uint64_t limit = (files == 1 ? 256 : 128) * mib;
        require(reader.get_sampled_plain_bytes() == files * limit, "sampling byte limit");
        require(kmers == files * (limit - 3 - k_length + 1), "truncated sequence counts");
    }
}

}
int main() {
    char pattern[] = "/tmp/fasta-reader-test-XXXXXX";
    auto path = ::mkdtemp(pattern);
    if (!path) return 1;
    int status = 0;
    try {
        for (auto k : {1U, 5U, 31U, 128U}) {
            k_length = k;
            for (auto block : {2048U, 4096U})
                for (auto threads : {1U, 3U}) check(path, threads, 4096, block);
        }
        check_sampling_limit(path);
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n'; status = 1;
    }
    std::filesystem::remove_all(path);
    return status;
}
