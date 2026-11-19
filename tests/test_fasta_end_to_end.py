"""Check format selection and counts against an independent sequence oracle."""
import collections
import gzip
import pathlib
import random
import subprocess
import sys
import tempfile

from high_prefix_format import read_high
from low_msd_format import read_low


def decode(record, k):
    return ''.join('ACGT'[(record[i // 4] >> (6 - 2 * (i % 4))) & 3] for i in range(k))


def expected_counts(sequences, k):
    result = collections.Counter()
    for sequence in sequences:
        sequence = sequence.upper()
        for i in range(len(sequence) - k + 1):
            word = sequence[i:i+k]
            if set(word) <= set('ACGT'):
                result[min(word[::-1], word.translate(str.maketrans('ACGT', 'TGCA')))] += 1
    return result


def write_input(path, sequences):
    fasta = path.name.removesuffix('.gz').endswith(('.fa', '.fasta'))
    if fasta:
        content = ''.join(f'>r{i}\r\n' + '\r\n'.join(seq[j:j+37] for j in range(0, len(seq), 37)) + '\r\n'
                          for i, seq in enumerate(sequences)).rstrip('\r\n')
    else:
        content = ''.join(f'@r{i}\n{seq}\n+\n' + '?' * len(seq) + '\n' for i, seq in enumerate(sequences))
    if path.suffix == '.gz':
        with gzip.open(path, 'wt') as out:
            out.write(content)
    else:
        path.write_text(content)


def run(binary, source, output, k, bloom):
    output.mkdir()
    return subprocess.run([binary, str(source), str(k), '10', '1', str(output),
                           '1024', '1', '4294967295', '65535', f'--bloom={bloom}'],
                          text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=60)


def main():
    binary = sys.argv[1]
    rng = random.Random(73)
    sequences = [''.join(rng.choices('ACGT', k=230)) for _ in range(12)]
    sequences += sequences[:4] * 3 + ['acgt' * 60, 'ACGT' * 40 + 'NNN' + 'TGCA' * 40, 'ACG']
    with tempfile.TemporaryDirectory(prefix='fasta-e2e-') as tmp:
        root = pathlib.Path(tmp)
        cases = []
        for suffix in ('fastq', 'fq', 'fastq.gz', 'fq.gz', 'fa', 'fasta', 'fa.gz', 'fasta.gz'):
            path = root / ('reads.' + suffix)
            write_input(path, sequences)
            cases.append((str(path), sequences))
        for suffixes in [('fq', 'fastq.gz', 'fq.gz'), ('fa', 'fasta.gz', 'fa.gz')]:
            directory = root / suffixes[0]
            directory.mkdir()
            parts = []
            for i, suffix in enumerate(suffixes):
                part = directory / f'{i}.{suffix}'
                write_input(part, sequences[i::3])
                parts.append(str(part))
            cases += [(','.join(parts), sequences), (str(directory), sequences)]
        for k in (5, 33, 128):
            expected = expected_counts(sequences, k)
            for mode in ('on', 'off'):
                for index, (source, _) in enumerate(cases):
                    output = root / f'out-{k}-{mode}-{index}'
                    result = run(binary, source, output, k, mode)
                    assert result.returncode == 0, result.stdout
                    actual = collections.Counter(decode(record, k) for record in read_low(output / 'low.bin', k))
                    for record, count in read_high(output / 'high.bin', k, 2):
                        actual[decode(record, k)] += count
                    assert actual == expected, (source, k, mode, result.stdout)
        empty = root / 'empty'
        empty.mkdir()
        for i, source in enumerate((str(empty), '', str(root / 'reads.fa') + ',' + str(root / 'reads.fq'), str(root / 'reads.txt'))):
            result = run(binary, source, root / f'bad-{i}', 5, 'off')
            assert result.returncode != 0, result.stdout
            assert any(message in result.stdout for message in ('No input files', 'Mixed FASTQ/FASTA', 'Unsupported input extension')), result.stdout
            assert not (root / f'bad-{i}' / 'infos.bin').exists(), result.stdout
    print('FASTA/FASTQ aliases, mixed compression, both pipelines and independent counts: passed')


if __name__ == '__main__':
    main()
