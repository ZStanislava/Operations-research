import argparse
import copy
from scipy.optimize import linprog
import traceback
from math import floor, ceil
import random as rand
import time

import warnings
warnings.filterwarnings("ignore")

rand.seed(11)

def input(filename):#считывание данных из файла
    f = open(filename, "r")
    type = f.readline()

    # размерность задачи 
    n = int(f.readline())
    # коэффицикнты в функционале
    c = list(map(int, f.readline().split()))

    A = []
    b = []
    with_no_additional = []

    # обработка ограничений типа равенство из файла 
    count = int(f.readline())
    for _ in range(count):
        line = list(map(int, f.readline().split()))
        A.append(line[:-1])
        b.append(line[-1])
        # избавляемся от отрицательности справа
        if b[-1] < 0:
            b[-1] *= -1
            A[-1] = list(map(lambda x: -x, A[-1]))
        with_no_additional.append(len(A) - 1)#индекс используется для отслеживания переменных, которые могут принимать дробные значения.

    shift = 0

    #  обработка ограничений типа больше из файла 
    count = int(f.readline())
    for _ in range(count):
        line = list(map(int, f.readline().split()))
        A.append(line[:-1])
        b.append(line[-1])
        # выравниваем на добавочные из прошлых условий
        for _ in range(shift):
            A[-1].append(0)
        # добавочная переменная
        A[-1].append(-1)
        shift += 1
        # выравниваем количество переменных
        c.append(0)
        for a in A[:-1]:
            a.append(0)
        # избавляемся от отрицательности справа
        if b[-1] < 0:
            b[-1] *= -1
            A[-1] = list(map(lambda x: -x, A[-1]))
        else:
            with_no_additional.append(len(A) - 1)

    # обработка ограничений типа меньше из файла 
    count = int(f.readline())
    for _ in range(count):
        line = list(map(int, f.readline().split()))
        A.append(line[:-1])
        b.append(line[-1])
        # выравниваем на добавочные из прошлых условий
        for _ in range(shift):
            A[-1].append(0)
        # добавочная переменная

        A[-1].append(1)
        shift += 1
        # выравниваем количество переменных
        c.append(0)
        for a in A[:-1]:
            a.append(0)
        # избавляемся от отрицательности справа
        if b[-1] < 0:
            b[-1] *= -1
            A[-1] = list(map(lambda x: -x, A[-1]))
            with_no_additional.append(len(A) - 1)

    # обработка количество отрицательных иксов из файла
    count = int(f.readline())
    # их индексы с нуля
    negative = list(map(int, f.readline().split())) if count > 0 else []
    for n in negative:
        c[n] *= -1
        for a in A:
            a[n] *= -1

    # обработка количество беззнаковых иксов-положительных  из файла 
    count = int(f.readline()) 
    # их индексы с нуля
    indefinite = list(map(int, f.readline().split())) if count > 0 else []
    for i in range(len(indefinite)):
        c.insert(indefinite[i] + i + 1, -c[indefinite[i] + i])
        for a in A:
            a.insert(indefinite[i] + i + 1, -a[indefinite[i] + i])

    # обработка количество целых иксов из файла
    count = int(f.readline())
    # их индексы с нуля
    integer = list(map(int, f.readline().split())) if count > 0 else []

    return A, b, c, with_no_additional, negative, indefinite, integer, n, type

def print_matrix(A):#принимает матрицу A и выводит ее на экран
    for line in A:
        for elem in line:
            print("{:10.4f}".format(elem), '\t', end='')
        print()
def print_matrix_and_vector(A, b):#принимает матрицу A и вектор b, затем выводит их на экран
    for line, b_elem in zip(A, b):
        for elem in line:
            print("{:10.4f}".format(elem), '\t', end='')
        print(" | {:10.4f}".format(b_elem))

def sub_line(fr, what, mult):#выполняется операция вычитания строки, умноженной на число, из другой строки.
    for i in range(len(fr)):
        fr[i] -= what[i] * mult

def div_line(line, divider):#выполняется операция деления строки на число.
    for i in range(len(line)):
        line[i] /= divider

def order_basic(A, b, basic):# текущая строка с базисной переменной перемещается в начало базисных строк.
    for index in range(len(basic)):
        for line in range(len(A)):
            if A[line][basic[index]] == 1:
                A[index], A[line] = A[line], A[index]
                b[index], b[line] = b[line], b[index]
                continue
            
def found_in_index(z_string, basic):# поиска индекса переменной, входящей в базис, с наименьшим коэффициентом в строке целевой функции.
    res = None
    variants = []
    for ind in range(len(z_string)):
        if ind in basic:
            continue
        if res is None or z_string[ind] < z_string[res]:
            res = ind
    return res

def step(z_string, A, b, basic, z_result):#один шаг симплекс-метода 
    out_string = None # ищет индекс строки, соответствующей выбору выходящей переменной (столбца) 
    if min(b) < 0: #выбираем вводимый столбец
        out_string = b.index(min(b))#Если в строке функции (строке Z) есть отрицательные коэффициенты, выбирается тот столбец, коэффициент в котором отрицателен и является наибольшим по модулю.
        in_index = None
        for ind in range(len(A[0])):
            
            if A[out_string][ind] >= 0 or ind in basic:
                continue
            if (in_index is None) or (abs(z_string[ind]/A[out_string][ind]) < abs(z_string[in_index]/A[out_string][in_index])):
                in_index = ind
    else:
        in_index = found_in_index(z_string, basic)
        for string in range(len(A)):#Если существуют положительные значения в строке функции (z_string), выбирается переменная с наименьшим коэффициентом в отношении к соответствующим коэффициентам в векторе b.
            if (A[string][in_index] > 0) and ((out_string is None) or ((b[string] / A[string][in_index]) < (b[out_string] / A[out_string][in_index]))):
                out_string = string

    assert out_string is not None, "extr not found"# выведены, если симплекс-метод не может найти подходящий столбец для ввода 
    assert in_index is not None, "no solution" #не может найти решение 
    
    basic.pop(out_string)
    basic.append(in_index)
    basic.sort()

    b[out_string] /= A[out_string][in_index]#Обновление выходной строки:Делит элементы выходной строки на значение элемента входящего столбца в выходной строке,
    div_line(A[out_string], A[out_string][in_index])#Делит всю выходную строку на значение элемента входящего столбца в выходной строке.

    z_result[0] -= b[out_string] * z_string[in_index]#Вычитает из текущего значения целевой функции произведение значения выходной строки и соответствующего коэффициента в строке целевой функции.
    sub_line(z_string, A[out_string], z_string[in_index])
    
    for string in range(len(A)):
        if string == out_string:
            continue
        
        b[string] -= b[out_string] * A[string][in_index]#Вычитает из соответствующего элемента вектора правых частей произведение значения выходной строки и соответствующего коэффициента в строке матрицы.
        sub_line(A[string], A[out_string], A[string][in_index])

    order_basic(A, b, basic)#упорядочивает базисные переменные, перемещая строки с базисными переменными в начало матрицы A и вектора b

def process_answers(x, negative, indefinite):#обработку и представление ответов после выполнения симплекс-метода.
    for n in negative:
        x[n] *= -1
    for ind in indefinite:
        first = x.pop(ind)
        second = x.pop(ind + 1)
        x.insert(ind, first - second)
    return x
    
def answers(A, b, basic, negative, indefinite):#обработку и представление ответов после выполнения симплекс-метода.
    x = [0] * len(A[0])
    for index in range(len(basic)):
        x[basic[index]] = b[index]
    return process_answers(x, negative, indefinite)
 
def simplex(A, b, c, z_result, negative, indefinite, n, type, basic=None, no_copy=False):
    if not no_copy:
        A = copy.deepcopy(A)
        b = b.copy()
        c = c.copy()
    
    negative = negative.copy()
    indefinite = indefinite.copy()

    if type == "min\n":
        c = list(map(lambda x: -x, c))
        
    z_string = list(map(lambda x: -x, c))

    if basic is None:
        basic = list(range(len(A[0]) - len(A), len(A[0])))
    
    while min(z_string) < 0 or min(b) < 0: #пока есть отрицательные значения в векторе целевой функции (z_string) или в векторе правых частей (b).
        step(z_string, A, b, basic, z_result)#На каждом шаге цикла вызывается функция step, которая реализует один шаг симплекс-метода.

    for i in range(len(z_string)):#После выхода из цикла осуществляется проверка, есть ли альтернативные решения 
        if i not in basic and z_string[i] == 0:#если есть переменные с нулевыми коэффициентами в строке целевой функции и они не входят в базис
            break#Если такая переменная с нулевым коэффициентом найдена, цикл завершается
        
    x = answers(A, b, basic, negative, indefinite)[0:n]
    
    if type == "min\n":
        z_result[0] *= -1# Если исходная задача на минимум, инвертируем знак значения целевой функции.
        
    return z_result[0], x, basic, z_string

def M_method(A, b, c, with_no_additional, negative, indefinite, n, type, no_copy=False):
    if not no_copy:
        A = copy.deepcopy(A)
        b = b.copy()
        c = c.copy()

    M = 1_000#значение M (большое положительное число), которое используется для создания искусственных переменных.

    if type == "min\n": #Если тип задачи "min", инвертируются знаки коэффициентов вектора c и добавляется искусственная переменная с отрицательным коэффициентом.
        c = list(map(lambda x: -x, c))
    
    for add in with_no_additional:#Это цикл по всем индексам переменные, которые могут принимать дробные значения).
        A[add].append(1)#В этой строке в матрицу A добавляется новая колонка
        c.append(-M)# вектор коэффициентов целевой функции c добавляется новый элемент, который соответствует коэффициенту искусственной переменной
        for string in range(len(A)):# цикл проходит по строкам матрицы A
            if string != add:#для каждой строки, кроме той, которую мы обозначили переменной add, 
                A[string].append(0)#добавляется новый элемент с значением 0 в конец этой строки

    z_string = list(map(lambda x: -x, c))

    z_result = [0]

    for add in with_no_additional:
        sub_line(z_string, A[add], M)
        z_result[0] -= b[add] * M # вычитания из вектора z_string произведения строки матрицы A (которая соответствует переменной, способной принимать дробные значения) на число M
        
    c = list(map(lambda x: -x, z_string))

    order_basic(A, b, range(len(A[0]) - len(A), len(A[0])))
    z, x, basic, z_string = simplex(A, b, c, z_result, negative, indefinite, n, "max\n", no_copy=no_copy)
    if type == "min\n":
        z *= -1
    return z, x, basic, z_string


def int_check(x, integer):#проверку на целочисленность
    res = []
    for ind in integer:
        if x[ind] - floor(x[ind]) > 10**(-4):
            res.append(ind)
    if len(res) == 0:
        return None
    return res[rand.randint(0, len(res) - 1)]


def round_array(x):#округление значений переменных.
    for i in range(len(x)):
        if x[i] - floor(x[i]) < 10**(-4):
            x[i] = floor(x[i])
        if ceil(x[i]) - x[i] < 10**(-4):
            x[i] = ceil(x[i])    


def gamory(A, b, c, with_no_additional, integer, negative, n, type):
    A = copy.deepcopy(A)
    b = b.copy()
    c = c.copy()
    count = 0
    if type == "min\n":
        c = list(map(lambda x: -x, c))
        
    if len(with_no_additional) == 0:#если нет переменных, которые могут принимать дробные значения (with_no_additional), вызывается функция simplex.
        z, x, basic, z_string = simplex(A, b, c, [0], negative, indefinite, n, "max\n", no_copy=True)
    else:#Иначе вызывается функция M_method для использования метода искусственного базиса.
        z, x, basic, z_string = M_method(A, b, c, with_no_additional, negative, indefinite, n, "max\n", no_copy=True)

    
    not_int_ind = int_check(x, integer)#Вызывается функция int_check для проверки целочисленности значений переменных. Если есть переменные с нецелочисленными значениями, начинается процесс модификации задачи.
    while not_int_ind is not None:#пока есть хотя бы одна нецелочисленная переменная
        count += 1#используется для отслеживания количества итераций цикла.
        
            
        line_number = basic.index(not_int_ind)#выбирается строка матрицы ограничений A, которая соответствует выбранной нецелочисленной переменной:
        b_frac = b[line_number] - floor(b[line_number])#Получаем дробную часть, вычитая целую часть из исходного значения
        new_line = []#новая строка new_line, которая используется для добавления новой искусственной переменной и обновления матрицы A и вектора b
        for ind in range(len(A[line_number])):#проходит по всем элементам строки матрицы ограничений A
            if ind in basic:# является ли текущий индекс ind базисным 
                new_line.append(0)#Если да, то добавляется 0 в новую строку.
            elif A[line_number][ind] > 0:#Если элемент A положителен, добавляется отрицательное значение в новую строку
                new_line.append(-A[line_number][ind])#если отрицателен, добавляется значение того элемента в новую строку. Это часть процесса модификации задачи для симплекс-метода.
            else:
                new_line.append((b_frac) / (1 - b_frac) * A[line_number][ind])#значение
        for a in A:
            a.append(0)
        z_string.append(0)
        new_line.append(1)
        A.append(new_line)
        b.append(-b_frac)
        
        basic.append(len(A[0]) - 1)
        order_basic(A, b, basic) #обновления базиса
        z_result = [z]
        
        c = map(lambda x: -x, z_string)
        z, x, basic, z_string = simplex(A, b, c, z_result, negative, indefinite, n, "max\n", no_copy=True, basic=basic)#вызывается симплекс-метод 
        
        if count % 100 == 0:
            print(z)
       
        round_array(x)
        round_array(b)
        not_int_ind = int_check(x, integer)

    if type == "min\n":
        z *= -1
    return z, x


def round_matrix(A):
    for i in range(len(A)):
        for j in range(len(A[i])):
            A[i][j] = round(A[i][j], 6)#округляется до 6 знаков после запятой с использованием функции round.
        

def calc_extr(A, b, c, with_no_additional, negative, indefinite, integer, n, type):
    print(type, end='')

    integrality=[1 if x in integer else 0 for x in range(0, len(A[0]))]
    if sum(integrality) == 0:
        integrality = None
        
    print("linprog:")

    if type == "min\n":
        result = linprog(c, A_eq=A, b_eq=b, integrality=integrality)
    else:
        result = linprog(list(map(lambda x: -x, c)), A_eq=A, b_eq=b, integrality=integrality)
        

    print("sucess: ", result.success)
    if result.success:
        print('{:.6f}'.format(result.fun if type == "min\n" else -result.fun))
        for x in process_answers(result.x, negative, indefinite)[0:n]:
            print('{:.6f}'.format(x), end=' ')
        print()
    print()

    try:
        if len(integer) > 0:
            print("Gamory")
            print(time.time())
            z, x = gamory(A, b, c, with_no_additional, integer, negative, n, type)
            print(time.time())
            print('{:.6f}'.format(z))
            for elem in x:
                print('{:.6f}'.format(elem), end=' ')
            print()
        elif len(with_no_additional) == 0:#сли нет переменных, которые могут принимать дробные значения (len(with_no_additional) == 0), выполняется расчет с использованием симплекс-метода (simplex).
            print("simplex")
            z, x, _, _ = simplex(A, b, c, [0], negative, indefinite, n, type)
            print('{:.6f}'.format(z))
            for elem in x:
                print('{:.6f}'.format(elem), end=' ')
            print()
        else:#В противном случае, выполняется расчет с использованием метода искусственного базиса (M_method).
            print("M-method")
            z, x, _, _ = M_method(A, b, c, with_no_additional, negative, indefinite, n, type)
            print('{:.6f}'.format(z))
            for elem in x:
                print('{:.6f}'.format(elem), end=' ')
            print()
    except Exception as e:#Обработка исключений:
        print(e)
    

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()
A, b, c, with_no_additional, negative, indefinite, integer, n, type = input(args.filename) 

if type == "extr\n":
    calc_extr(A, b, c, with_no_additional, negative, indefinite, integer, n, "min\n")
    print()
    calc_extr(A, b, c, with_no_additional, negative, indefinite, integer, n, "max\n")
else:
    calc_extr(A, b, c, with_no_additional, negative, indefinite, integer, n, type)
    
    

    

    
