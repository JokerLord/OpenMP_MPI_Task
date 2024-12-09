# Реализация задания по курсу "Суперкомпьютерное моделирование и технологии" 2 курса магистратуры ВМК МГУ

Сборка и запуск на IBM Polus (MPI версия):
```bash
mpixlC -o main.out main.cpp -std=c++11 -O3
mpisubmit.pl -p <кол-во MPI процессов> main.out --stdout std.out --stderr std.err 128 20 1. 1. 1. 0.01
```
