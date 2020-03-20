import random
import numpy as np

x1_min = 20
x1_max = 70
x2_min = 5
x2_max = 40
x3_min = 20
x3_max = 45
xm_min = (x1_min + x2_min + x3_min) / 3
xm_max = (x1_max + x2_max + x3_max) / 3
y_min = 200 + xm_min
y_max = 200 + xm_max

kohrenTable = {1: 0.9065, 2: 0.7679, 3: 0.6841, 4: 0.6287, 5: 0.5892,
               6: 0.5598, 7: 0.5365, 8: 0.5175, 9: 0.5017, 10: 0.4884}

studentTable = {4: 2.776, 8: 2.306, 12: 2.179, 16: 2.120, 20: 2.086, 24: 2.064, 28: 2.048}

fisherTable = {
    1: {1: 164.4, 2: 199.5, 3: 215.7, 4: 224.6, 5: 230.2, 6: 234.0, 12: 244.9, 24: 249.0},
    2: {1: 18.5, 2: 19.2, 3: 19.2, 4: 19.3, 5: 19.3, 6: 19.3, 12: 19.4, 24: 19.4},
    3: {1: 10.1, 2: 9.6, 3: 9.3, 4: 9.1, 5: 9.0, 6: 8.9, 12: 8.7, 24:  8.6},
    4: {1: 7.7, 2: 6.9, 3: 6.6, 4: 6.4, 5: 6.3, 6: 6.2, 12: 5.9, 24: 5.8},
    5: {1: 6.6, 2: 5.8, 3: 5.4, 4: 5.2, 5: 5.1, 6: 5.0, 12: 4.7, 24: 4.5},
    6: {1: 6.0, 2: 5.1, 3: 4.8, 4: 4.5, 5: 4.4, 6: 4.3, 12: 4.0, 24: 3.8},
    7: {1: 5.5, 2: 4.7, 3: 4.4, 4: 4.1, 5: 4.0, 6: 3.9, 12: 3.6, 24: 3.4},
    8: {1: 5.3, 2: 4.5, 3: 4.1, 4: 3.8, 5: 3.7, 6: 3.6, 12: 3.3, 24: 3.1},
    12:{1: 4.8, 2: 3.9, 3: 3.5, 4: 3.3, 5: 3.1, 6: 3.0, 12: 2.7, 24: 2.5}}

x_norm = [[-1, -1, -1],
          [-1, 1, 1],
          [1, -1, 1],
          [1, 1, -1]]

x_arr = [[20, 5, 20],
         [20, 40, 45],
         [70, 5, 45],
         [70, 40, 20]]

m = 2
y_arr = [[random.randint(int(y_min), int(y_max)) for i in range(m)] for j in range(4)]

print("X1   X2  X3  Y1   Y2")
for i in range(4):
    print(x_arr[i], y_arr[i])


def kohren(dispersion, m, gt):
    gp = max(dispersion) / sum(dispersion)
    return gp < gt[m - 1]


def student(dispersion_reproduction, m, y_aver, xn):
    dispersion_statistic_mark = (dispersion_reproduction / (4 * m)) ** 0.5

    beta = [1 / 4 * sum(y_aver[j] for j in range(4))]
    for i in range(3):
        b = 0
        for j in range(4):
            b += y_aver[j] * xn[j][i]
        beta.append(1 / 4 * b)

    t = []
    for i in beta:
        t.append(abs(i) / dispersion_statistic_mark)

    return t[0] > studentTable[(m - 1) * 4], t[1] > studentTable[(m - 1) * 4], \
           t[2] > studentTable[(m - 1) * 4], t[3] > studentTable[(m - 1) * 4]


def fisher(m, d, y_aver, yo, dispersion_reproduction, fisherTable):
    dispersion_ad = 0
    for i in range(4):
        dispersion_ad += (yo[i] - y_aver[i]) ** 2

    dispersion_ad = dispersion_ad * m / (4 - d)

    fp = dispersion_ad / dispersion_reproduction

    return fp < fisherTable[(m - 1) * 4][4 - d]


"""У внутрішньому while (92строка) перевіряється однорідність дисперсій,
а якщо вони неоднорідні => m=m+1 та знову розраховуються У та дисперсії, що залежать від m,
цикл працює доки дисперсії не стануть однорідними"""

"""У зовнішньому while (91строка) перевіряється адекватність математичної моделі оригіналу.
Якщо модель неадекватна, або всі коефіцієнти статистично незначущі за Стьюдентом => m+=1
цикл працює доки адекватність не підтвердиться критерієм Фішера"""
while True:
    while True:
        y_aver = []
        for i in range(4):
            y_aver.append(sum(y_arr[i]) / m)

        dispersion = []
        for i in range(len(y_arr)):
            dispersion.append(0)
            for j in range(m):
                dispersion[i] += (y_aver[i] - y_arr[i][j]) ** 2
            dispersion[i] /= m

        dispersion_reproduction = sum(dispersion) / 4

        if kohren(dispersion, m, kohrenTable):
            break
        else:
            m += 1
            for i in range(4):
                y_arr[i].append(random.randint(int(y_min), int(y_max)))

    k = student(dispersion_reproduction, m, y_aver, x_norm)
    d = sum(k)

    mx1 = (x_arr[0][0] + x_arr[1][0] + x_arr[2][0] + x_arr[3][0]) / 4
    mx2 = (x_arr[0][1] + x_arr[1][1] + x_arr[2][1] + x_arr[3][1]) / 4
    mx3 = (x_arr[0][2] + x_arr[1][2] + x_arr[2][2] + x_arr[3][2]) / 4
    my = sum(y_aver) / 4

    a11 = (x_arr[0][0] ** 2 + x_arr[1][0] ** 2 + x_arr[2][0] ** 2 + x_arr[3][0] ** 2) / 4
    a22 = (x_arr[0][1] ** 2 + x_arr[1][1] ** 2 + x_arr[2][1] ** 2 + x_arr[3][1] ** 2) / 4
    a33 = (x_arr[0][2] ** 2 + x_arr[1][2] ** 2 + x_arr[2][2] ** 2 + x_arr[3][2] ** 2) / 4
    a12 = (x_arr[0][0] * x_arr[0][1] + x_arr[1][0] * x_arr[1][1] + x_arr[2][0] * x_arr[2][1] + x_arr[3][0] * x_arr[3][1]) / 4
    a13 = (x_arr[0][0] * x_arr[0][2] + x_arr[1][0] * x_arr[1][2] + x_arr[2][0] * x_arr[2][2] + x_arr[3][0] * x_arr[3][2]) / 4
    a23 = (x_arr[0][1] * x_arr[0][2] + x_arr[1][1] * x_arr[1][2] + x_arr[2][1] * x_arr[2][2] + x_arr[3][1] * x_arr[3][2]) / 4

    a1 = (x_arr[0][0] * y_aver[0] + x_arr[1][0] * y_aver[1] + x_arr[2][0] * y_aver[2] + x_arr[3][0] * y_aver[3]) / 4
    a2 = (x_arr[0][1] * y_aver[0] + x_arr[1][1] * y_aver[1] + x_arr[2][1] * y_aver[2] + x_arr[3][1] * y_aver[3]) / 4
    a3 = (x_arr[0][2] * y_aver[0] + x_arr[1][2] * y_aver[1] + x_arr[2][2] * y_aver[2] + x_arr[3][2] * y_aver[3]) / 4

    b_arr = []
    b_arr.append(np.linalg.det(np.array([[my, mx1, mx2, mx3],
                                         [a1, a11, a12, a13],
                                         [a2, a12, a22, a23],
                                         [a3, a13, a23, a33]]))/np.linalg.det(np.array([[1, mx1, mx2, mx3],
                                                                                        [mx1, a11, a12, a13],
                                                                                        [mx2, a12, a22, a23],
                                                                                        [mx3, a13, a23, a33]])))
    b_arr.append(np.linalg.det(np.array([[1, my, mx2, mx3],
                                         [mx1, a1, a12, a13],
                                         [mx2, a2, a22, a23],
                                         [mx3, a3, a23, a33]])) / np.linalg.det(np.array([[1, mx1, mx2, mx3],
                                                                                          [mx1, a11, a12, a13],
                                                                                          [mx2, a12, a22, a23],
                                                                                          [mx3, a13, a23, a33]])))
    b_arr.append(np.linalg.det(np.array([[1, mx1, my, mx3],
                                         [mx1, a11, a1, a13],
                                         [mx2, a12, a2, a23],
                                         [mx3, a13, a3, a33]])) / np.linalg.det(np.array([[1, mx1, mx2, mx3],
                                                                                          [mx1, a11, a12, a13],
                                                                                          [mx2, a12, a22, a23],
                                                                                          [mx3, a13, a23, a33]])))
    b_arr.append(np.linalg.det(np.array([[1, mx1, mx2, my],
                                         [mx1, a11, a12, a1],
                                         [mx2, a12, a22, a2],
                                         [mx3, a13, a23, a3]])) / np.linalg.det(np.array([[1, mx1, mx2, mx3],
                                                                                          [mx1, a11, a12, a13],
                                                                                          [mx2, a12, a22, a23],
                                                                                          [mx3, a13, a23, a33]])))
    """b_arr = ([b_arr[i] * k[i] for i in range(4)])"""
    """Попередній рядок можна роздокументовати, це відкине незначущі коефіцієнти, але збільшить похибку"""

    yo = []
    for i in range(4):
        yo.append(b_arr[0] + b_arr[1] * x_arr[i][0] + b_arr[2] * x_arr[i][1] + b_arr[3] * x_arr[i][2])

    if d == 4:
        m += 1
        for i in range(4):
            y_arr[i].append(random.randint(int(y_min), int(y_max)))

    elif fisher(m, d, y_aver, yo, dispersion_reproduction, fisherTable):
        break
    else:
        m += 1
        for i in range(4):
            y_arr[i].append(random.randint(int(y_min), int(y_max)))

'''checking'''
checks = []
errors = 0
for i in range(4):
    print("Y average", i, " = ", y_aver[i])
    print("Y found  ", i, " = ", yo[i])
    checks.append(round(y_aver[i], 10) == round(yo[i], 10))
for j in range(4):
    if checks[j] == False:
        errors += 1
        print("Our test failed!")

if errors == 0:
    print("Successful check!")
