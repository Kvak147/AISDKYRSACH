import random 
import time
import matplotlib.pyplot as plt
import numpy as np
# Функции для работы с матрицами

Size = 14

def Determinate(A):
    det = 0
    n = len(A[0])
    if (n == 2):
        return (A[0][0]*A[1][1] - A[0][1]*A[1][0])
    i = 0
    for i in range(n):
        b = []
        for j in range(n):
            if j == 0:
                continue
            else:
                string = []
                for k in range(n):
                    if k != i:
                        string.append(A[j][k])
            b.append(string)
        if (1+i+1)%2 == 0:
            det = det + A[0][i]*Determinate(b)
        else:
            det = det + (-1)*A[0][i]*Determinate(b)
    return det

def Transpose(A):
    n = len(A)
    B = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            B[j][i] = A[i][j]
    return B

def TransposeV(A):
    n = len(A)
    B = [[0] for _ in range(n)]
    for i in range(n):
        B[i][0] = A[i]
    return B

def TransposeC(A):
    n = len(A)
    B = [0 for _ in range(n)]
    for i in range(n):
        B[i] = A[i][0]
    return B

def Reverse(A):
    n = len(A)
    det = Determinate(A)
    if det == 0:
        print("Singular matrix")
        return 
    A = Transpose(A)
    for i in range(n):
        for j in range(n):
            A[i][j] = A[i][j]/det
    return A

def MultiplicationM(A,B):
    if len(A[0]) != len(B):
        print("Error len(A[0]) != len(B)")
        return
    else:
        C = [[0 for _ in range(len(B[0]))] for _ in range(len(A))]
        for i in range(len(A)):
            for j in range(len(B[0])):
                s = 0
                for k in range(len(A[0])):
                    s += A[i][k]*B[k][j]
                C[i][j] = s
        return C

def MultiplicationVectors(A,B):
    C = [[0 for _ in range(len(B))] for _ in range(len(B))]
    for i in range(len(B)):
        for j in range(len(B)):
            C[i][j] = A[i][0]*B[j]
    return C

def AdditionM(A,B):
    if len(A) != len(B) or len(A[0]) != len(B[0]):
        print("len(A) != len(B) or len(A[0]) != len(B[0]")
        return
    else:
        n = len(A)
        m = len(A[0])
        for i in range(n):
            for j in range(m):
                A[i][j] += B[i][j]
        return A

def MultiplicationV(A,B):
    n = len(A)
    for i in range(n):
        for j in range(n):
            A[i][j] *= B
    return A

def AdditionV(A,B):
    n = len(A)
    for i in range(n):
        for j in range(n):
            A[i][j] += B
    return A

def Triangle(A):
    n = len(A)
    for i in range(1,n):
        for j in range(0,i):
            k1 = A[j][j]
            k2 = A[i][j]
            for z in range(n):
                A[i][z] = A[i][z]*k1 - A[j][z]*k2
    return A

def Sign(A):
    if A > 0:
        return 1
    elif A < 0:
        return -1
    else:
        return 0

def Norma(A):
    n = len(A)
    s = 0
    for i in range(n):
        s += A[i]**2
    return s**(0.5)

def NormaColumn(A):
    n = len(A)
    s = 0
    for i in range(n):
        s += A[i][0]**2
    return s**(0.5)

def SubtractV(A,B):
    n = len(A)
    for i in range(n):
        A[i] = A[i] - B[i]
    return A

def VectorValueMultip(A,B):
    n = len(A)
    for i in range(n):
        A[i] *= B
    return A

def ColumnValueMultip(A,B):
    n = len(A)
    for i in range(n):
        A[i][0] *= B
    return A

def StrColumnMultiplucation(A,B):
    s = 0
    n = len(A)
    for i in range(n):
        s += A[i]*B[i][0]
    return s

def MethodHouseHolder(A,V):
    n = len(A)
    R = A.copy()
    Q = [[0 if i!=j else 1 for i in range(n)] for j in range(n)]
    for k in range(n):
        x = [R[i+k][k] if i+k < n else 0 for i in range(n-k)]
        x = TransposeV(x)
        sign = Sign(x[0][0])
        normx = NormaColumn(x)
        e = TransposeV([1 if i == 0 else 0 for i in range(n-k)])
        e[0][0] = normx * sign
        u = AdditionM(x,e)
        u = ColumnValueMultip(u,1/NormaColumn(u))
        H = [[0 if i!=j else 1 for i in range(n)] for j in range(n)]
        um = MultiplicationVectors(u,TransposeC(u))
        M = [[0 for _ in range(n)] for _ in range(n)]
        i = 0
        while (k+i < n):
            j = 0
            while (k+j < n):
                M[i+k][j+k] = um[i][j]
                j += 1
            i+=1
        H = AdditionM(H, MultiplicationV(M,-2))
        R = MultiplicationM(H,R)
        Q = MultiplicationM(Q,Transpose(H))
        V = MultiplicationM(V,Q)
    return Q,R,V

def CheckSxodimost(A):
    s = 0
    n = len(A)
    for i in range(n):
        for j in range(i+1,n):
            s += A[i][j]**2
    s = s**(0.5)
    if s <= 10**(-10):
        return True
    else:
        return False

def QR(A):
    k = 10000
    n = len(A)
    V = [[0 if i!=j else 1 for i in range(n)] for j in range(n)]
    for i in range(k):
        Q, R, V = MethodHouseHolder(A,V)
        A = MultiplicationM(R,Q)
        if CheckSxodimost(A):
            print(i)
            break
    return A, V
# ------------------------------
def GenerationD(n):
    A = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        A[0][i] = float(float(random.randint(-15,15)) + float(random.randint(-15,15))/100)
        while A[0][i] == 0:
            A[0][i] = float(float(random.randint(-15,15)) + float(random.randint(-15,15))/100)
    for i in range(1,n):
        for j in range(n):
            if j%2 == 0:
                A[i][j] = A[i-1][j]+i/10
            else:
                A[i][j] = A[i-1][j]-i/10
            while A[i][j] == 0:
                A[i][j] = float(float(random.randint(-15,15)) + float(random.randint(-15,15))/100)
    return A

def GenerationR(n):
    A = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            A[i][j] = float(float(random.randint(-15,15)) + float(random.randint(-15,15))/100)
            while A[i][j] == 0:
                A[i][j] = float(float(random.randint(-15,15)) + float(random.randint(-15,15))/100)
    return A

def ShowVectors(V):
    n = len(V)
    print("Собственные вектора:")
    for i in range(n):
        print("Вектор",i+1)
        for j in range(n):
            print(V[j][i])
    return


def ShowValueOfMatrix(A):
    n = len(A)
    print("Собственные числа:")
    for i in range(n):
        print("λ",i+1,") - ",A[i][i])
    return 

def ShowMatrix(A):
    n = len(A)
    print("Матрица:")
    for i in range(n):
        for j in range(n):
            print(A[i][j],"\t", end='')
        print()

def ShowQRMatrixs(cn,cm):
    for i in range (2,cm):
        if cn == 1:
            A = GenerationD(i)
        else:
            A = GenerationR(i)
        ShowMatrix(A)
        Start = time.time()
        A, V = QR(A)
        End = time.time()
        ShowMatrix(A)
        ShowValueOfMatrix(A)
        ShowVectors(V)
        print("Время выполнения - ", End - Start)

def ShowQRforMatrix():
    print("Какого размера матрица?")
    n = int(input())
    A = [[0 for _ in range(n)] for _ in range(n)]
    print("Вводите по очереди значение каждого элемента, идя построчно")
    for i in range(n):
        for j in range(n):
            A[i][j] = float(input())
    ShowMatrix(A)
    Start = time.time()
    A, V = QR(A)
    End = time.time()
    ShowMatrix(A)
    ShowValueOfMatrix(A)
    ShowVectors(V)
    print("Время выполнения - ", End - Start)

def DrawGraphQRRandom():
    times = [0 for i in range(Size)]
    SizesMatrixs =[0 for i in range(Size)]
    for i in range(Size):
        SizesMatrixs[i] = i+2
        A = GenerationR(i+2)
        Start = time.time()
        A, V = QR(A)
        End = time.time()
        times[i] = End - Start
    coeff = np.polyfit(SizesMatrixs, times,3)
    poly = np.poly1d(coeff)
    Sizes_new = np.linspace(min(SizesMatrixs),max(SizesMatrixs),100)
    times_new = poly(Sizes_new)
    print(coeff)
    plt.scatter(SizesMatrixs,times,label="Time of выполнение",color='purple')
    plt.plot(Sizes_new, times_new, color='red', label='Cubic Regression')
    plt.xlabel('Size of matrix')
    plt.ylabel('Time')
    plt.title('Время выполнения QR алгориthm на случайной матрице')
    plt.legend()
    plt.grid(True)
    plt.xticks(np.arange(min(SizesMatrixs), max(SizesMatrixs) + 1, 1.0))
    plt.savefig('graphRandom.png')
    plt.show()
    return 

def DrawGraphQRDeterminated():
    times = [0 for i in range(Size)]
    SizesMatrixs =[0 for i in range(Size)]
    for i in range(Size):
        SizesMatrixs[i] = i+2
        A = GenerationD(i+2)
        Start = time.time()
        A, V = QR(A)
        End = time.time()
        times[i] = End - Start
    coeff = np.polyfit(SizesMatrixs, times,3)
    poly = np.poly1d(coeff)
    Sizes_new = np.linspace(min(SizesMatrixs),max(SizesMatrixs),100)
    times_new = poly(Sizes_new)
    print(coeff)
    plt.scatter(SizesMatrixs,times,label="Time of выполнение",color='purple')
    plt.plot(Sizes_new, times_new, color='red', label='Cubic Regression')
    plt.xlabel('Size of matrix')
    plt.ylabel('Time')
    plt.title('Время выполнения QR алгориthm на почти вырожденной матрице')
    plt.legend()
    plt.grid(True)
    plt.xticks(np.arange(min(SizesMatrixs), max(SizesMatrixs) + 1, 1.0))
    plt.savefig('graphDeterminated.png')
    plt.show()
    return 

def Main():
    x = 1
    while x != 0:
        print("Что вы хотите?")
        print("0 - выйти\n1 - сделать QR алгоритм для случайных матриц до определенного размера n\n2 - сделать QR алгоритм для почти вырожденных матриц до определенного размера n\n3 - Ввести свою матрицу для QR алгоритма\n4 - вывести график для случайных матриц\n5 - вывести график для почти вырожденных матриц")
        x = int(input())
        if x == 1:
            print("До какого размера идём?")
            n = int(input())
            ShowQRMatrixs(x,n+1)
        elif x == 2:
            print("До какого размера идём?")
            n = int(input())
            ShowQRMatrixs(x,n+1)
        elif x == 3:
            ShowQRforMatrix()
        elif x == 4:
            DrawGraphQRRandom()
        elif x ==  5:
            DrawGraphQRDeterminated()
        else:
            print("Ввели что-то не то)))")
    return

Main()