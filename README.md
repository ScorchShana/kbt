# kbt

`kbt` counts canonical k-mers in FASTQ files. A k-mer is a sequence of `k`
bases. A sequence and its reverse complement share one count.

The program is written in C++20 for Linux and is designed for large datasets
and machines with large amounts of RAM. It reads plain or gzip-compressed FASTQ
files and saves compact binary results. Two tools turn those results into text:

- `dump_tool` exports k-mer sequences and their counts.
- `histogram_tool` reports how many distinct k-mers have each count.

## Build

You need GCC with C++20 support, CMake 3.14 or newer, and Make or Ninja.
NUMA support is enabled when the libnuma headers and library are available.
For example, on Ubuntu or Debian:

```sh
sudo apt install build-essential cmake libnuma-dev
```

From the project directory, initialize the zlib-ng submodule and build the three
programs:

```sh
git submodule update --init --recursive
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target kbt dump_tool histogram_tool --parallel 4
```

The executables are `build/kbt`, `build/dump_tool`, and `build/histogram_tool`.
The build uses `-march=native` to optimize for the build machine's CPU.

## Quick start

Create a directory for the binary results, then count 31-mers using a thread
setting of 16 and a memory budget of 64 GiB:

```sh
mkdir -p results
./build/kbt reads.fastq.gz 31 16 64 results
```

By default, `kbt` keeps k-mers with a count of at least 2 and caps each count at
255. Export the saved results as a sorted, tab-separated file:

```sh
./build/dump_tool precise results 4 1 kmers.tsv 1 4294967295 --sort
```

Create a histogram of the saved counts:

```sh
./build/histogram_tool results histogram.tsv 1 255
```

To keep k-mers that occur only once as well, use a separate result directory
and set `filter_min` to 1:

```sh
mkdir -p results-all
./build/kbt reads.fastq.gz 31 16 64 results-all 1024 1
./build/dump_tool precise results-all 4 1 kmers-all.tsv 1 4294967295 --sort
./build/histogram_tool results-all histogram-all.tsv 1 255
```

## Counting options

```text
kbt <fastq_file> <k_length> <n_thread> <memory_limit_gb> <res_dir>
    [map_capacity] [filter_min] [filter_max] [count_max] [mapping_filename]
    [--bloom=on|off] [--first-cache=auto|off] [--first-reclaim=spill]
```

The command above is wrapped for readability. Supply the positional arguments
in the order shown. To set a later optional argument, supply the earlier ones
too.

| Required argument | Meaning |
| --- | --- |
| `fastq_file` | One FASTQ file, a comma-separated list of files, or an input directory. Plain FASTQ and `.gz` files are supported. |
| `k_length` | Number of bases per k-mer, from 1 to 128. |
| `n_thread` | Thread setting used to divide work between processing stages. Must be at least 6. |
| `memory_limit_gb` | Memory budget as a whole number of GiB. One GiB is 1,073,741,824 bytes. Thread stacks and other process overhead also use memory. |
| `res_dir` | An existing directory for working files and final binary results. |

| Optional argument | Default | Meaning |
| --- | --- | --- |
| `map_capacity` | `1024` | Maximum capacity setting for per-node hash tables. Rounded up to a power of two, with a minimum of 1024. |
| `filter_min` | `2` | Lowest stored count to keep, inclusive. Use 1 to include singletons. |
| `filter_max` | `4294967295` | Highest stored count to keep, inclusive. |
| `count_max` | `255` | Count cap, from 1 to 4,294,967,295. Higher counts are stored as this value. |
| `mapping_filename` | `mapping.tmp` | Path used for temporary disk-backed memory. Relative paths use the current working directory. |

Counts are capped **before** frequency filtering. For example, with
`count_max=255`, a k-mer seen 300 times has a stored count of 255. It is kept
when the filter range includes 255.

To keep counts from 2 through 100 and raise the cap to 65,535:

```sh
mkdir -p results-range
./build/kbt reads.fastq.gz 31 16 64 results-range 1024 2 100 65535
```

To retain singletons and use the full 32-bit count range:

```sh
mkdir -p results-wide
./build/kbt reads.fastq.gz 31 16 64 results-wide 1024 1 4294967295 4294967295
```

### Input files

Multiple files contribute to the same set of counts. For example, paired-end
files can be counted together:

```sh
./build/kbt reads_R1.fastq.gz,reads_R2.fastq.gz 31 16 64 results
```

You can also pass a directory:

```sh
./build/kbt fastq-inputs/ 31 16 64 results
```

Directory input includes every non-hidden entry, so keep only input FASTQ
files in that directory. Use standard four-line FASTQ records. Uppercase and
lowercase bases are accepted, and `U` is treated as `T`. A base such as `N`
breaks a k-mer window; counting resumes after enough valid bases have been read.
Windows do not cross read boundaries.

### Processing modes

| Option | Default | Meaning |
| --- | --- | --- |
| `--bloom=on` or `--bloom=off` | `on` | Use the Bloom-based pipeline, or send every k-mer directly to the counting tree. Both produce the same retained counts. |
| `--first-cache=auto` or `--first-cache=off` | `auto` | In Bloom mode, enable or disable the in-memory cache for first occurrences. |
| `--first-reclaim=spill` | `spill` | Write cached first occurrences to temporary storage when memory needs to be reclaimed. |

For example, to count without the Bloom filter:

```sh
mkdir -p results-no-bloom
./build/kbt reads.fastq.gz 31 16 64 results-no-bloom --bloom=off
```

## Binary results

After a successful run, the result directory contains:

| File | Contents |
| --- | --- |
| `infos.bin` | Metadata: k-mer length and count cap. |
| `high.bin` | Retained k-mers from the counting tree, with their stored counts. This file can also contain singletons. |
| `low.bin` | Retained singletons from the Bloom pipeline. Each record has an implicit count of 1. |

The two k-mer files contain separate keys. Read both for the complete saved
result. With the default `filter_min=2`, `low.bin` is empty. With
`--bloom=off`, all retained k-mers are in `high.bin` and `low.bin` is empty.

Keep these three files together for export and histogram generation. Files
named `low.pending.bin.*` and `*.final.tmp` are working files used during
counting; a successful run removes or renames them. Use a separate result
directory for each run, since final files in that directory are replaced.

## Export k-mers as text

```text
dump_tool <precise|approximate> <res_dir> <max_threads> <max_memory_gb>
          <output_file> [min_freq=1] [max_freq=10000] [--sort]
```

- `precise` reads `high.bin` and includes `low.bin` when the requested range
  contains count 1. Use this for the complete saved result.
- `approximate` reads only `high.bin`, leaving out any singletons in `low.bin`.
- `max_threads` sets the number of export workers; `max_memory_gb` is the
  export memory budget in GiB and can be fractional.
- `min_freq` and `max_freq` select an inclusive range of stored counts.
- `--sort` writes sequences in alphabetical order. Without it, output order
  is unspecified. Sorting can use temporary disk space in the result directory.

For example, export every saved count with four workers and a 1 GiB budget:

```sh
./build/dump_tool precise results-all 4 1 kmers.tsv 1 4294967295 --sort
```

The file has two tab-separated columns and no header:

```text
AAAAA	3
ACGTA	1
CCCCC	2
```

Each row gives a canonical k-mer sequence and its stored count. In this example,
`AAAAA` and its reverse complement `TTTTT` occur three times in total.

The default export maximum is 10,000. Pass a larger `max_freq`, as above, when
your saved counts exceed that value. Export filters apply to the results kept
by the counting step.

## Create a frequency histogram

```text
histogram_tool <res_dir> <output_file> [min_freq=1] [max_freq=10000]
```

This tool reads the saved binary results and writes one row for each frequency
in the requested inclusive range, including frequencies with no k-mers. It has
no mode, thread-count, or memory-budget arguments.

```sh
./build/histogram_tool results-all histogram.tsv 1 255
```

The output has two tab-separated columns and no header. For the three k-mers
shown above, the first five rows would be:

```text
1	1
2	1
3	1
4	0
5	0
```

The first column is the stored frequency. The second is the number of distinct
canonical k-mers with that frequency. For example, `3  1` means one distinct
k-mer has count 3; it does not mean three distinct k-mers.

To include singletons in the histogram, keep them during counting with
`filter_min=1`. A count cap of 255 puts every k-mer seen at least 255 times in
the frequency-255 bin.

## Console output

`kbt` prints the input settings, the thread split, and a final summary. The main
summary fields are:

| Field | Meaning |
| --- | --- |
| `Total read k-mer count` | Number of valid k-mer occurrences read, including repeated occurrences and those later filtered out. |
| `High frequency unique k-mer` | Number of distinct retained keys written to `high.bin`. Add `Low frequency records kept` for the total number of distinct saved k-mers. |
| `Low frequency records scanned` | First-occurrence records read back from temporary storage. |
| `Low frequency records recovered` | Scanned records added back to the counting tree. |
| `FIRST drain singleton records` | Singletons found while draining the in-memory first-occurrence cache. |
| `Low frequency records kept` | Singletons written to `low.bin`. |
| `Low frequency records discarded` | Singletons omitted from the low-frequency output by the counting filter. |
| `Unique k-mers filtered below filter_min` | Number of distinct keys removed by the lower frequency bound. |
| `Unique k-mers filtered above filter_max` | Number of distinct keys removed by the upper frequency bound. |
| `Pending encoded bytes written` / `read` | Encoded temporary data written and read during counting. |
| `Total elapsed us` | Counting run time in microseconds. Divide by 1,000,000 for seconds. |

The low-frequency processing fields are zero when `--bloom=off` is used.
`dump_tool` and `histogram_tool` print progress separately from their result
files.

For additional timing, memory, and queue statistics, build with diagnostics:

```sh
cmake -S . -B build-diagnostics -DCMAKE_BUILD_TYPE=Release -DTEST_MODE=ON
cmake --build build-diagnostics --target kbt --parallel 4
```
