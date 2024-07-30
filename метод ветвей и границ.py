import argparse
import copy
from scipy.optimize import linprog
from math import ceil, floor

import warnings
warnings.filterwarnings("ignore")

best = -2**31
best_vector = None

class SimplexTable:# представления и работы с таблицей симплекс-метода 
    def __init__(self, A = [], b = [], c = []):
        self.A = copy.deepcopy(A)
        self.b = b.copy()
        self.c = c.copy()

    # return type of optimization
    def read(self, filename) -> str:
        self.basic = []

        f = open(filename, "r")
        type = f.readline()
    
        # размерность задачи 
        self.n = int(f.readline())
        # коэффицикнты в функционале
        self.c = list(map(int, f.readline().split()))

        self.A = []
        self.b = []
        self.with_no_additional = []

        # ограничений типа равенство
        count = int(f.readline())
        for _ in range(count):
            line = list(map(int, f.readline().split()))
            self.A.append(line[:-1])
            self.b.append(line[-1])
            # избавляемся от отрицательности справа
            if self.b[-1] < 0:
                self.b[-1] *= -1
                self.A[-1] = list(map(lambda x: -x, self.A[-1]))
            self.with_no_additional.append(len(self.A) - 1)

        shift = 0

        # ограничений типа больше
        count = int(f.readline())
        for _ in range(count):
            line = list(map(int, f.readline().split()))
            self.A.append(line[:-1])
            self.b.append(line[-1])
            # выравниваем на добавочные из прошлых условий
            for _ in range(shift):
                self.A[-1].append(0)
            # добавочная переменная
            self.A[-1].append(-1)
            shift += 1
            # выравниваем количество переменных
            self.c.append(0)
            for a in self.A[:-1]:
                a.append(0)
            # избавляемся от отрицательности справа
            if self.b[-1] < 0:
                self.b[-1] *= -1
                self.A[-1] = list(map(lambda x: -x, self.A[-1]))
                self.basic.append(len(self.A[-1]) - 1)
            else:
                self.with_no_additional.append(len(self.A) - 1)

        # ограничений типа меньше
        count = int(f.readline())
        for _ in range(count):
            line = list(map(int, f.readline().split()))
            self.A.append(line[:-1])
            self.b.append(line[-1])
            # выравниваем на добавочные из прошлых условий
            for _ in range(shift):
                self.A[-1].append(0)
            # добавочная переменная

            self.A[-1].append(1)
            shift += 1
            # выравниваем количество переменных
            self.c.append(0)
            for a in self.A[:-1]:
                a.append(0)
            self.basic.append(len(self.A[-1]) - 1)
            # избавляемся от отрицательности справа
            if self.b[-1] < 0:
                self.b[-1] *= -1
                self.A[-1] = list(map(lambda x: -x, self.A[-1]))
                self.with_no_additional.append(len(self.A) - 1)

        # количество отрицательных иксов
        count = int(f.readline())
        # их индексы с нуля
        self.negative = list(map(int, f.readline().split())) if count > 0 else []
        for n in self.negative:
            self.c[n] *= -1
            for a in self.A:
                a[n] *= -1

        # количество беззнаковых иксов
        count = int(f.readline()) 
        # их индексы с нуля
        self.indefinite = list(map(int, f.readline().split())) if count > 0 else []
        for i in range(len(self.indefinite)):
            self.c.insert(self.indefinite[i] + i + 1, -self.c[self.indefinite[i] + i])
            for a in self.A:
                a.insert(self.indefinite[i] + i + 1, -a[self.indefinite[i] + i])

        # количество целых иксов
        count = int(f.readline())
        # их индексы с нуля
        self.integer = list(map(int, f.readline().split())) if count > 0 else []

        self.z_string = list(map(lambda x: -x, self.c))
        self.result = 0

        return type

    def print(self):
        for line, b_elem in zip(self.A, self.b):
            for elem in line:
                print("{:10.4f}".format(elem), '\t', end='')
            print(" | {:10.4f}".format(b_elem))

    def div_line(self, index, divider):
        for i in range(len(self.A[index])):
            self.A[index][i] /= divider
        self.b[index] /= divider

    def sub_line(self, reduced, subtracted, multiplier):
        for i in range(len(self.A[reduced])):
            self.A[reduced][i] -= self.A[subtracted][i] * multiplier
        self.b[reduced] -= self.b[subtracted] * multiplier

    def sub_z_line(self, subtracted, multiplier):
        for i in range(len(self.z_string)):
            self.z_string[i] -= self.A[subtracted][i] * multiplier
        self.result -= self.b[subtracted] * multiplier

    def set_basic(self, basic = None):#устанавливает базис переменных таблицы
        if basic is None:
            self.basic = list(range(len(self.A[0]) - len(self.A), len(self.A[0])))
        else:
            self.basic = basic.copy()
    
    def order_basic(self):
        for index in range(len(self.basic)):
            for line in range(len(self.A)):
                if self.A[line][self.basic[index]] == 1:
                    self.A[index], self.A[line] = self.A[line], self.A[index]
                    self.b[index], self.b[line] = self.b[line], self.b[index]
                    continue

    def update_basic(self, out_string, in_index):#обновляет базис после выполнения одного шага симплекс-метод
        self.basic.pop(out_string)
        self.basic.append(in_index)
        self.basic.sort()

        self.div_line(out_string, self.A[out_string][in_index])
        for string in range(len(self.A)):
            if string != out_string:
                self.sub_line(string, out_string, self.A[string][in_index])
        self.sub_z_line(out_string, self.z_string[in_index])

        self.order_basic()
        
    def step(self):#шаг симплекс-метода. 
        if not(self.is_valid() and not self.is_optimum()):#Проверяется, что текущее решение допустимо (is_valid()) и не является оптимальным (not is_optimum()).
            print("solution should be valid but not optimal")
        assert self.is_valid() and not self.is_optimum(), "solution should be valid but not optimal"
        in_index = self.z_string.index(min(self.z_string))#Выбор переменной для ввода в базис 
        if not(in_index not in self.basic):#Проверяется, что переменная с минимальным коэффициентом в строке целевой функции (z_string) не входит уже в базис
            print("it is forbidden to introduce into the basis what is already in the basis")
        assert in_index not in self.basic, "it is forbidden to introduce into the basis what is already in the basis"
        
        out_string = None
        for string in range(len(self.A)):
            if self.A[string][in_index] > 0 and \
            (out_string is None or self.b[string] / self.A[string][in_index] < self.b[out_string] / self.A[out_string][in_index]):
                out_string = string

        assert out_string is not None, "extr not found"
        self.update_basic(out_string, in_index)

    def dual_step(self): 
        if not (self.is_optimum() and not self.is_valid()):#Проверяется, что текущее решение является оптимальным и не допустимым
            print("solution should be optional, but not valid")
        assert self.is_optimum() and not self.is_valid(), "solution should be optional, but not valid"#Если условие не выполняется, возникает исключение
        out_string = self.b.index(min(self.b))#Находится выводимая строка- ищем индекс строки с минимальным значением вектора правых частей b. 
        in_index = None
        for ind in range(len(self.A[out_string])):
            if self.A[out_string][ind] >= 0 or ind in self.basic:#Переменная уже входит в базис или ее коэффициент в строке out_string неотрицателен?
                continue#Проверка на минимальность отношения коэффициента переменной в строке целевой функции к коэффициенту переменной в строке out_string.
            if (in_index is None) or \
                (abs(self.z_string[ind] / self.A[out_string][ind]) < abs(self.z_string[in_index] / self.A[out_string][in_index])):
                
                in_index = ind

        assert in_index is not None, "no solution" #Проверяется, что переменная in_index не осталась равной None
        self.update_basic(out_string, in_index)#обновляет базис

    def is_optimum(self):#возвращает True, если минимальный элемент в векторе z_string (коэффициенты целевой функции) неотрицателен,
        return (min(self.z_string) >= 0)
    
    def is_valid(self):#Возвращает True, если минимальный элемент в векторе b (вектор правых частей) неотрицателен
        return (min(self.b) >= 0)

    def process_answers(self, x):#преобразует вектор x
        for n in self.negative:
            x[n] *= -1
        for ind in self.indefinite:
            first = x.pop(ind)
            second = x.pop(ind + 1)
            x.insert(ind, first - second)
        return x

    def get_answers(self):#ответы симплекс-метода на текущей итераци
        x = [0] * len(self.A[0])
        for index in range(len(self.basic)):
            x[self.basic[index]] = self.b[index]
        return self.process_answers(x)[0:self.n]

    def int_check(self):# выполняет проверку целочисленности переменных в текущем базисе. 
        x = self.get_answers()
        res = None
        diff = 0
        for ind in self.integer:
            if ind not in self.basic: 
                assert x[ind] == 0
                continue
            up = ceil(x[ind]) - x[ind]
            down = x[ind] - floor(x[ind])
            if min(up, down) < 10**(-6):
                continue
            if res is None or down > diff:
                res = self.basic.index(ind)
                diff = down
        return res, diff
    
    def get_first_float(self):#возвращает индекс первой переменной, которая не является целым числом.
        x = self.get_answers()
        for ind in self.integer:
            if x[ind] != int(x[ind]):#если значение не является целым числом, функция возвращает индекс соответствующей переменной в базисе (self.basic).
                return self.basic.index(ind)
            up = ceil(x[ind]) - x[ind]
            down = x[ind] - floor(x[ind])
            if min(up, down) > 10**(-6):
                return self.basic.index(ind)
    
    def drop_M_vars(self):#выполняет удаление переменных, добавленных при использовании метода искусственного базиса (M-метод).
        count = len(self.with_no_additional)
        for _ in range(count):
            self.z_string.pop()
            self.c.pop()
            for a in self.A:
                a.pop()

def simplex(table: SimplexTable, type = "max\n"):
    mult = 1
    if type == 'min\n':
        table = copy.deepcopy(table)
        table.c, table.z_string = table.z_string, table.c
        mult = -1
    
    while not table.is_optimum():
        table.step()
    return table.result * mult, table.get_answers()
    
def M_method(table: SimplexTable, type = "max\n"):
    M = 1_000_000

    mult = 1
    if type == 'min\n':
        table = copy.deepcopy(table)
        table.c, table.z_string = table.z_string, table.c
        mult = -1

    for add in table.with_no_additional:
        table.basic.append(len(table.c))
        table.A[add].append(1)
        table.c.append(-M)
        table.z_string.append(M)
        for string in range(len(table.A)):
            if string != add:
                table.A[string].append(0)
    
    for add in table.with_no_additional:
        table.sub_z_line(add, M)

    table.basic.sort()
    table.order_basic()

    result, x = simplex(table)

    return result * mult, x

def gomory(table: SimplexTable, type):
    mult = 1
    if type == 'min\n':
        table = copy.deepcopy(table)
        table.c, table.z_string = table.z_string, table.c
        mult = -1
    
    if len(table.with_no_additional) == 0:
        simplex(table)
    else:
        M_method(table)
        table.drop_M_vars()

    line_number, diff = table.int_check()
    while line_number is not None:
        new_line = []
        
        for ind in range(len(table.A[line_number])):
            if ind in table.basic:
                new_line.append(0)
            elif table.A[line_number][ind] > 0:
                new_line.append(-table.A[line_number][ind])
            else:
                new_line.append(diff / (1 - diff) * table.A[line_number][ind])
        for a in table.A:
            a.append(0)
        table.z_string.append(0)
        table.c.append(0)
        new_line.append(1)
        table.A.append(new_line)
        table.b.append(-diff)

        table.basic.append(len(table.A[line_number]) - 1)
        table.order_basic()

        assert not table.is_valid(), "after gomory step should be invalid"
        while not table.is_valid():
            table.dual_step()

        line_number, diff = table.int_check()      

    return table.result * mult, table.get_answers()
    
def border_step(table: SimplexTable):
    while not table.is_valid():#пока решение недопустимо 
        table.dual_step()#выполняем двойственный симплекс-метод
    
    global best, best_vector
    if table.result <= best: #Если текущий результат в table меньше или равен глобальному лучшему результату (best), функция завершается
        return
    line = table.get_first_float()#Поиск первой переменной с плавающей точкой: Если все переменные являются целыми числами, line устанавливается в None.
    if line is None: #Проверка на целочисленное решение:Если нет переменных с плавающей точкой (line равно None)
        if table.result > best:#является ли текущий результат лучше глобального лучшего результата. 
            best = table.result#Если да, она обновляет глобальный лучший результат
            best_vector = table.get_answers()# и лучший вектор 

        return
    
    less = copy.deepcopy(table)#Функция создает копии less текущей таблицы table
    less_line = [-x for x in table.A[line]] + [1]# новой строке less_line коэффициенты переменных из строки line умножаются на -1 (отрицательная копия), а затем добавляется 1 в конец строки. 
    less_line[table.basic[line]] = 0#устанавливается значение элемента в столбце basic, соответствующего переменной в строке line, в 0.

    for a in less.A:
        a.append(0)#Для каждой строки a в матрице A текущей таблицы less, добавляется новый элемент (0) 
    less.A.append(less_line)# строка less_line добавляется в конец матрицы A
    
    less.basic.append(len(less_line) - 1)#добавляет индекс последней добавленной переменной (новой переменной, введенной в методе граничного поиска) в список basic.
    less.z_string.append(0)#добавляет коэффициент 0 для соответствующей новой переменной в вектор z_string.
    less.b.append(floor(table.b[line]) - table.b[line])#добавляет новое значение вектора правых частей less.b. как разность между нижней границей (округленной вниз) текущего значения table.b[line] и самим этим значением.

    try:#обработать исключение, возникающее при вызове рекурсивной функции border_step(less)

        border_step(less)
    except Exception as e: 

        pass

    more = copy.deepcopy(table)#Функция создает копии more текущей таблицы table
    more_line = [x for x in table.A[line]] + [1]
    more_line[table.basic[line]] = 0

    for a in more.A:
        a.append(0)
    more.A.append(more_line)
    
    more.basic.append(len(more_line) - 1)
    more.z_string.append(0)
    more.b.append(table.b[line] - ceil(table.b[line]))#Добавляет новое значение вектора правых частей, вычисленное как разность между текущим значением table.b[line] и его потолочным (округленным вверх) значением.

    try:

        border_step(more)
    except Exception as e:

        pass
    
def border(table: SimplexTable, type = "max\n"):
    mult = 1#фактором умножения для коррекции значения результата в зависимости от типа задачи (максимизация или минимизация).
    if type == 'min\n':
        table = copy.deepcopy(table)
        table.c, table.z_string = table.z_string, table.c#преобразовать задачу минимизации в эквивалентную задачу максимизации
        mult = -1#что результат будет умножен на -1. 

    if len(table.with_no_additional) == 0:#проверяет, имеются ли переменные добавочные в таблице. Если их нет, вызывается функция simplex(table)
        simplex(table)
    else:
        M_method(table)#иначе вызывается функция M_method(table)
        table.drop_M_vars()#затем удаляются добавочные переменные с помощью table.drop_M_vars().

    try:
        border_step(table)
    except Exception as e:

        pass

    assert best_vector is not None, "no solution" #лучшее решение (best_vector) не является None,
    return mult * best, best_vector  #и возвращает оптимальное значение и вектор переменных, умноженные на mult= 1 или -1 в зависимости от типа

def process_answers(x, negative, indefinite):#
    for n in negative:
            x[n] *= -1
    for ind in indefinite:
        first = x.pop(ind)
        second = x.pop(ind + 1)
        x.insert(ind, first - second)
    return x

def calc_extr(table: SimplexTable, type):
    global best, best_vector
    best = -2**31
    best_vector = None

    table = copy.deepcopy(table)
    print(type, end='')

    integrality=[1 if x in table.integer else 0 for x in range(0, len(table.A[0]))]
    if sum(integrality) == 0:
        integrality = None
    
    print("linprog:")
    if type == "min\n":
        result = linprog(table.c, A_eq=table.A, b_eq=table.b, integrality=integrality)
    else:
        result = linprog(list(map(lambda x: -x, table.c)), A_eq=table.A, b_eq=table.b, integrality=integrality)
    print("sucess: ", result.success)
    if result.success:
        print('{:.6f}'.format(result.fun if type == "min\n" else -result.fun))
        for x in process_answers(result.x, table.negative, table.indefinite)[0:table.n]:
            print('{:.6f}'.format(x), end=' ')
        print()

    try:
        if len(table.integer) > 0:
            print("Border")
            z, x = border(table, type)
            print('{:.6f}'.format(z))
            for elem in x:
                print('{:.6f}'.format(elem), end=' ')

        elif len(table.with_no_additional) == 0:
            print("simplex")
            z, x = simplex(table, type)
            print(z, x)
        else:
            print("M-method")
            z, x = M_method(table, type)
            print(z, x)
    except Exception as e:
        print(e)

parser = argparse.ArgumentParser()
parser.add_argument('filename')
table = SimplexTable()
filename ="D:/Новая папка/иатэ/4 курс/ИО/мои/2/11b.txt"
type = table.read(filename)

if type == "extr\n":
    calc_extr(table, "min\n")
    print()
    calc_extr(table, "max\n")
else:
    calc_extr(table, type)
    
    
