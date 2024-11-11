# mpi-inversed-matrix
>Блочное нахождение обратной матрицы с разделённой памятью методом Жордана без выделения дополнительной памяти с выбором главного элемента по всей матрице

Алгоритм с виртуализацией присоединенной матрицы, позволяющий обойтись одной матрицей описан тут https://file.scirp.org/pdf/AM_2013100413422038.pdf

## Установка mpi на Ubuntu
>sudo apt-get -y install mpi-default-dev

Компиляция: mpicxx -isystem /opt/impi-5.1.3.223/intel64/include <флаги> <файлы>
## Флаги компиляции
-pg -O3 -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format
## Параметры запуска
Параметры командной строки:
>mpirun -np p ./a.out n m r s filename
0) p - количество процессов
1) n – размерность матрицы,
2) m – размер блока,
3) r – количество выводимых значений в матрице,
4) s – задает номер формулы для инициализации матрицы, 0 при вводе
матрицы из файла
5) filename – имя файла, откуда надо прочитать матрицу. Этот аргумент отсутствует,
если s!= 0.

Пример запуска:
>mpirun -np 4 ./a.out 2000 90 10 3

Запуск на 4 mpi-процессах матрицы размера 2000x2000, использовать блочный алгоритм с размером блока 90, выводить не более 10 строк и столбцов, заполнять матрицу по формуле 3
## Пример запуска
>![output](https://github.com/user-attachments/assets/285ef15b-c7e5-4d08-8c54-c3f469119cb3)

T1 - время нахождения обратной матрицы
T2 - время нахождения невязок
Невязки (при n < 11000): ![equation](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D%5Cbg%7Bwhite%7D$r_%7B1%7D=%7C%7CAA%5E%7B-1%7D-E%7C%7C_%7B1%7D$,$r_%7B2%7D=%7C%7CA%5E%7B-1%7DA-E%7C%7C_%7B1%7D$)

При n = 18000, m = 60, ускорение на 24 процессах к 1 процессу ~23.3
