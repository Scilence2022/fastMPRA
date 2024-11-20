# fastMPRA

fastMPRA is a high-performance tool for analyzing Massively Parallel Reporter Assay (MPRA) data. It processes FASTQ files to identify and analyze specific sequence regions of interest.

## Features

- Fast processing of single-end and paired-end FASTQ files
- Multi-threaded implementation for high performance
- Configurable k-mer size for sequence matching
- Support for multiple design sequences and zones
- Flexible configuration through INI files
- Quality-aware read merging for paired-end data

## Installation

### Prerequisites

- GCC compiler
- zlib development files
- pthread support

### Building

```bash
git clone https://github.com/Scilence2022/fastMPRA.git
cd fastMPRA
make

## Usage

bash
./fastMPRA [-h] [-c config.ini] [-1 <read1.fq>] [-2 <read2.fq>] [-k kmer_size] [-t threads]

### Options

- `-h`: Show help message
- `-c <file>`: Config file (default: config.ini)
- `-1 <file>`: Input FASTQ file for read 1 (required)
- `-2 <file>`: Input FASTQ file for read 2 (optional, for paired-end data)
- `-k <int>`: K-mer size (default: 31)
- `-t <int>`: Number of threads (default: 5)
```

## Configuration File Format

The program uses an INI-style configuration file. Here's an example:

```ini
; FastMPRA-UH Configuration File
; 这是注释


[General]
k = 21                 
threads = 5               
min_overlap = 12          
max_mismatch_ratio = 0.1  

[Input]
fq1 = reads_1.fq.gz
fq2 = reads_2.fq.gz

[Output]
output_file = output.txt  
output_format = 0         

[Advanced]
batch_size = 10000        
quality_threshold = 20   


[Design]
label = DNA-design
seq = ccaggggtccccaataattacgatttaaatttgacataataatacgactcactataggGCTATGGTTAGTTCCCACGTTccagctcccatgtaggcgtgcccaaacACCTTAGTAGGTACTACTACAACGTCGCACCGACTACGTTATAACGGGACGCCACAGAGACTTTGTTAAGGCCCGTGGTAAGAACAATTACCGATTACCCACCCCTTTACTGCCATGTGTGACATTAGCAAGAGTCCAATCCCCCCGAAAGCTTAGTCGGTTAGTCCCAGGCACTTC

[ZONE]
name = T7-pro
start = 46
end = 50

[ZONE]
name = UMI
start = 80
end = 106

```

## Output Format

The program outputs tab-separated values with the following columns:
1. Design label
2. Zone name
3. Extracted sequence
4. Left k-mer found (+/-) 
5. Right k-mer found (+/-)

## Performance

The tool uses multiple optimization strategies:
- Multi-threading for parallel processing
- Efficient k-mer matching algorithms
- Batch processing of reads
- Memory-efficient data structures

## Authors

- Lifu Song [songlf@tib.cas.cn]
- Rongxing Wang
- Mei Wang

## Contact

For questions and support, please contact:
- Email: [songlf@tib.cas.cn](mailto:songlf@tib.cas.cn)

## License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.

## Citation

If you use fastMPRA in your research, please cite:
[Citation information to be added after publication]

