# hse22_hw1 - Седов Сергей 202


## I. Работаем в терминале

### 1. Заходим на сервер по ключу.
```console
ssh sasedov@92.242.58.92 -p 5222 -i github-manjaro-ssh
```

### 2. Создаём ссылку на папку с файлами в своей директории.
```console
ln -s /usr/share/data-minor-bioinf/assembly
```

### 3. Случайно выбираем чтения для paired-end и mate-pairs
```console
seqtk sample -s119 ./assembly/oil_R1.fastq 5000000 > sub_oil_R1.fastq
seqtk sample -s110 ./assembly/oil_R2.fastq 5000000 > sub_oil_R2.fastq
seqtk sample -s110 ./assembly/oilMP_S4_L001_R1_001.fastq 1500000 > sub_oilMP_S4_L001_R1_001.fastq
seqtk sample -s110 ./assembly/oilMP_S4_L001_R2_001.fastq 1500000 > sub_oilMP_S4_L001_R2_001.fastq
mkdir fastq
mv *.fastq fastq
```

### 4. Оцениваем качество исходных чтений с помощью FastQC
```console
mkdir fastqc
ls fastq/* | xargs -P 4 -tI{} fastqc -o fastqc {}
```

### 5. Общая статистика по файлам .fastqc через MultiQC
```console
mkdir multiqc
multiqc -o multiqc fastqc
```

### 6. Подрезаем чтения и удаляем адаптеры, переместим новые файлы в отдельную директорию
```console
platanus_trim sub_oil_*
platanus_internal_trim sub_oilMP*
mkdir fastq_trimmed
mv -v *trimmed fastq_trimmed
```

### 7. Запускаем FastQC и MultiQC для них
```console
mkdir fastqc_trimmed
ls fastq_trimmed/* | xargs -P 4 -tI{} fastqc -o fastqc_trimmed {}
mkdir multiqc_trimmed
multiqc -o multiqc_trimmed fastqc_trimmed
```

### 8. Собираем контиги из подреззаных чтений
```console
time platanus assemble -o Poil -t 8 -m 28 -f fastq_trimmed/sub_oil_R1.fastq.trimmed fastq_trimmed/sub_oil_R2.fastq.trimmed 2> contig.log
```

### 9. Теперь собираем скаффолды
```console
time platanus scaffold -o Poil -t 8 -c Poil_contig.fa -IP1 fastq_trimmed/sub_oil_R1.fastq.trimmed fastq_trimmed/sub_oil_R2.fastq.trimmed -OP2 fastq_trimmed/sub_oilMP_S4_L001_R1_001.fastq.int_trimmed fastq_trimmed/sub_oilMP_S4_L001_R2_001.fastq.int_trimmed 2> scaffold.log
```

### 10. Уменьшаем количество гэпов
```console
time platanus gap_close -o Poil -t 8 -c Poil_scaffold.fa -IP1 fastq_trimmed/sub_oil_R1.fastq.trimmed fastq_trimmed/sub_oil_R2.fastq.trimmed -OP2 fastq_trimmed/sub_oilMP_S4_L001_R1_001.fastq.int_trimmed fastq_trimmed/sub_oilMP_S4_L001_R2_001.fastq.int_trimmed 2> gapclose.log
```

### 11. Удаляем ненужные файлы с сервера
```console
rm -rf fastq
rm -rf fastq_trimmed
```

### 12. Копируем отчёты и файлы от platanus'а с сервера на локальную машину
```console
scp -i github-manjaro-ssh -r -P 5222 sasedov@92.242.58.92:/home/sasedov/multiqc/ /home/sergey/bioinf-minor/ 
scp -i github-manjaro-ssh -r -P 5222 sasedov@92.242.58.92:/home/sasedov/multiqc_trimmed/ /home/sergey/bioinf-minor/ 
scp -i github-manjaro-ssh -P 5222 sasedov@92.242.58.92:/home/sasedov/Poil_contig.fa /home/sergey/bioinf-minor/ 
scp -i github-manjaro-ssh -P 5222 sasedov@92.242.58.92:/home/sasedov/Poil_scaffold.fa /home/sergey/bioinf-minor/ 
scp -i github-manjaro-ssh -P 5222 sasedov@92.242.58.92:/home/sasedov/Poil_gapClosed.fa /home/sergey/bioinf-minor/ 
```

## Посмотрим отчёты MultiQC
### Для первоначальных чтений
![](https://github.com/NoMoreActimel/hse22_hw1/blob/main/data/multiqc_1.png)
![](https://github.com/NoMoreActimel/hse22_hw1/blob/main/data/multiqc_2.png)
![](https://github.com/NoMoreActimel/hse22_hw1/blob/main/data/multiqc_3.png)

### Для обрезанных чтений
![](https://github.com/NoMoreActimel/hse22_hw1/blob/main/data/multiqc_trimmed_1.png)
![](https://github.com/NoMoreActimel/hse22_hw1/blob/main/data/multiqc_trimmed_2.png)
![](https://github.com/NoMoreActimel/hse22_hw1/blob/main/data/multiqc_trimmed_3.png)

## II. Переходим в ipynb на ![Google Colab](https://colab.research.google.com/drive/11fBYrsHDWA9uT8SF-US3_YCFGUOw_5QM?usp=sharing)

### 1. Открываем файлы
```python
poil_contig_file = open('Poil_contig.fa', 'r')
poil_scaffold_file = open('Poil_scaffold.fa', 'r')
poil_gapClosed_file = open('Poil_gapClosed.fa', 'r')
```

### 2. Основная функция
```python
def get_file_info(file):
    """Basic anylysis for fasta file"""
    """:param file: input fasta file"""
    """:return: (number_of_sequences, total_length, max_length, longest_seq_start, N50)"""

    number_of_seq = 0
    cur_seq_start = 0
    cur_seq_length = 0
    total_length_read = 0
    sequences: list[tuple[int, int]] = []
    n50_length = 0

    for line in file.readlines():
        if line[0] == '>':
            if cur_seq_length:
                sequences.append((cur_seq_length, cur_seq_start))
            cur_seq_start = total_length_read
            cur_seq_length = 0
        else:
            cur_seq_length += len(line.strip())
        total_length_read += len(line)  
    
    if total_length_read:
        sequences.append((cur_seq_length, cur_seq_start))
    
    sequences.sort(reverse=True)

    max_length = sequences[0][0]
    longest_seq_start = sequences[0][1]

    number_of_seq = len(sequences)
    total_length = sum(sequence[0] for sequence in sequences)

    top_n_sum_length = 0
    for length, start in sequences:
        top_n_sum_length += length
        if top_n_sum_length >= total_length / 2:
            n50_length = length
            break
    
    return number_of_seq, total_length, max_length, longest_seq_start, n50_length
```

### 3. Пишем обёртку для вывода и получаем информацию для контигов и скаффолдов
```python
def print_file_info(file, file_info):
    print(f"Info about {file.name}:")
    print(f"\tNumber of sequences: {file_info[0]}")
    print(f"\tTotal length of sequences: {file_info[1]}")
    print(f"\tMax sequence length: {file_info[2]}")
    print(f"\tN50: {file_info[4]}\n")
```
```python
contig_info = get_file_info(poil_contig_file)
print_file_info(poil_contig_file, contig_info)
```
```text
Info about Poil_contig.fa:
	Number of sequences: 628
	Total length of sequences: 3926785
	Max sequence length: 179307
	N50: 52802
```

```python
scaffold_info = get_file_info(poil_scaffold_file)
print_file_info(poil_scaffold_file, scaffold_info)
```
```text
Info about Poil_scaffold.fa:
	Number of sequences: 73
	Total length of sequences: 3876446
	Max sequence length: 3834891
	N50: 3834891
```

### 4. Функция для подсчёта гэпов в самой длинной последовательности
```python
def count_gaps_for_max_seq(file, longest_seq_start, longest_seq_length):
    """
    Counts gaps and overall gaps' length in the longest sequence 
    :param file: input fasta file
    :param longest_seq_start: index where the longest sequence starts
    :param longest_seq_length: length of longest sequence
    :return: number_of_gaps, gaps_total_length
    """
    gap_cnt = 0
    cnt_N = 0
    cnt_N_in_a_row = 0

    file.seek(longest_seq_start)
    line = file.readline()
    assert line[0] == '>'

    sequence = file.read(longest_seq_length)
    for i, c in enumerate(sequence.strip()):
        if c == '\n':
            continue
        if c == 'N':
            cnt_N_in_a_row += 1
            cnt_N += 1
        elif cnt_N_in_a_row:
            gap_cnt += 1
            cnt_N_in_a_row = 0
    
    return gap_cnt, cnt_N
```

### 5. Вызываем для скаффолдов с обычными гэпами и сокращёнными
```python
scaffold_gap_info = count_gaps_for_max_seq(
    file=poil_scaffold_file, 
    longest_seq_start=scaffold_info[3], 
    longest_seq_length=scaffold_info[2]
)
print(f"Gaps for biggest scaffold: {scaffold_gap_info[0]}")
print(f"Their total length: {scaffold_gap_info[1]}")
```
```text
Gaps for biggest scaffold: 58
Their total length: 6430
```

```python
gapClosed_info = get_file_info(poil_gapClosed_file)
poil_gapClosed_file.seek(0)
gapClosed_gap_info = count_gaps_for_max_seq(
    file=poil_gapClosed_file, 
    longest_seq_start=gapClosed_info[3],
    longest_seq_length=gapClosed_info[2]
)
print(f"Gaps for biggest scaffold with gap_close: {gapClosed_gap_info[0]}")
print(f"Their total length: {gapClosed_gap_info[1]}")
```
```text
Gaps for biggest scaffold with gap_close: 8
Their total length: 1423
```

### 6. Создаём longest.fa
```python
with open('longest.fa', 'w') as longest_file:
    poil_scaffold_file.seek(scaffold_info[3])
    sequence = poil_scaffold_file.readline()
    assert sequence[0] == '>'
    sequence += poil_scaffold_file.read(scaffold_info[2])
    longest_file.write(sequence)    
```
