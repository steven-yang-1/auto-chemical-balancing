import numpy
import math
from numpy.linalg import matrix_rank
from scipy.linalg import null_space

MOL_MASS = {
    'H': 1.008,
    'He': 4.003,
    'Li': 6.941,
    'Be': 9.012,
    'B': 10.811,
    'C': 12.011,
    'N': 14.007,
    'O': 15.999,
    'F': 18.999,
    'Ne': 20.183,
    'Na': 22.990,
    'Mg': 24.305,
    'Al': 26.982,
    'Si': 28.086,
    'P': 30.974,
    'S': 32.07,
    'Cl': 35.453,
    'Ar': 39.948,
    'K': 39.092,
    'Ca': 40.078,
    'Sc': 44.956,
    'Ti': 47.867,
    'V': 50.942,
    'Cr': 51.996,
    'Mn': 54.938,
    'Fe': 55.845,
    'Co': 58.933,
    'Ni': 58.693,
    'Cu': 63.546,
    'Zn': 65.409,
    'Ga': 69.723,
    'Ge': 72.61,
    'As': 74.922,
    'Se': 78.96,
    'Br': 79.904,
    'Kr': 83.798,
    'Rb': 85.468,
    'Sr': 87.621,
    'Y': 88.906,
    'Zr': 91.224,
    'Nb': 92.907
}

TOKEN_ELEMENT = 0
TOKEN_NUMBER = 1
TOKEN_PLUS = 2
TOKEN_SPLIT = 3

STATUS_ALPHA = 0
STATUS_DIGIT = 1
STATUS_PLUS = 2
STATUS_EMPTY = 3
STATUS_SPLIT = 4


def lexer(in_str):
    buffer = []
    tokens = []
    status = STATUS_ALPHA
    for i, ch in enumerate(in_str):
        if status == STATUS_ALPHA:
            if ch.isalpha():
                if i > 0 and not "".join([buffer[0], ch]) in MOL_MASS:
                    if not "".join(buffer) in MOL_MASS:
                        raise Exception("化学反应式中出现了非法的元素，请仔细检查。")
                    tokens.append((TOKEN_ELEMENT, "".join(buffer)))
                    buffer = [ch]
                else:
                    buffer.append(ch)
            elif ch.isdigit():
                if not "".join(buffer) in MOL_MASS:
                    raise Exception("化学反应式中出现了非法的元素，请仔细检查。")
                tokens.append((TOKEN_ELEMENT, "".join(buffer)))
                buffer = [ch]
                status = STATUS_DIGIT
            elif ch == ' ' or ch == '\t':
                tokens.append((TOKEN_ELEMENT, "".join(buffer)))
                status = STATUS_EMPTY
            elif ch == '-' or ch == '>':
                status = STATUS_SPLIT
            elif ch == '+':
                tokens.append((TOKEN_ELEMENT, "".join(buffer)))
                status = STATUS_PLUS
            else:
                raise Exception("语法错误")
        elif status == STATUS_DIGIT:
            if ch.isalpha():
                tokens.append((TOKEN_NUMBER, int("".join(buffer))))
                buffer = [ch]
                status = STATUS_ALPHA
            elif ch.isdigit():
                buffer.append(ch)
            elif ch == ' ' or ch == '\t':
                tokens.append((TOKEN_NUMBER, int("".join(buffer))))
                status = STATUS_EMPTY
            elif ch == '-' or ch == '>':
                status = STATUS_SPLIT
            elif ch == '+':
                tokens.append((TOKEN_NUMBER, int("".join(buffer))))
                status = STATUS_PLUS
            else:
                raise Exception("语法错误")
        elif status == STATUS_EMPTY:
            if ch.isalpha():
                buffer = [ch]
                status = STATUS_ALPHA
            elif ch.isdigit():
                buffer = [ch]
                status = STATUS_DIGIT
            elif ch == '-' or ch == '>':
                status = STATUS_SPLIT
            elif ch == '+':
                status = STATUS_PLUS
        elif status == STATUS_PLUS:
            if ch.isalpha():
                buffer = [ch]
                status = STATUS_ALPHA
            elif ch.isdigit():
                buffer = [ch]
                status = STATUS_DIGIT
            elif ch == ' ' or ch == '\t':
                status = STATUS_EMPTY
            else:
                raise Exception("语法错误")
            tokens.append((TOKEN_PLUS, "+"))
        elif status == STATUS_SPLIT:
            if ch.isalpha():
                buffer = [ch]
                tokens.append((TOKEN_SPLIT, "->"))
                status = STATUS_ALPHA
            elif ch.isdigit():
                tokens.append((TOKEN_SPLIT, "->"))
                buffer = [ch]
                status = STATUS_DIGIT
            elif ch == ' ' or ch == '\t':
                tokens.append((TOKEN_SPLIT, "->"))
                buffer = []
                status = STATUS_EMPTY

    if len(buffer) > 0:
        if status == STATUS_ALPHA:
            if not "".join(buffer) in MOL_MASS:
                raise Exception("化学反应式中出现了非法的元素，请仔细检查。")
            tokens.append((TOKEN_ELEMENT, "".join(buffer)))
        elif status == STATUS_DIGIT:
            tokens.append((TOKEN_NUMBER, int("".join(buffer))))
        else:
            raise Exception("语法错误")

    return tokens


def extract_unique(tokens):
    exists = {}
    for current_token in tokens:
        if current_token[0] == TOKEN_ELEMENT:
            exists[current_token[1]] = True
    keys = exists.keys()
    return sorted(keys)


def format_result(tokens, coefs):
    str_buf = []
    cursor = 0
    for i in range(len(coefs)):
        coefs[i] = round(coefs[i])
    if coefs[cursor] != 1:
        str_buf.append(str(round(coefs[cursor])))
    cursor = cursor + 1
    for current_token in tokens:
        if current_token[0] == TOKEN_ELEMENT:
            str_buf.append(current_token[1])
        elif current_token[0] == TOKEN_NUMBER:
            str_buf.append(str(current_token[1]))
        elif current_token[0] == TOKEN_PLUS:
            str_buf.append(" + ")
            if coefs[cursor] != 1:
                str_buf.append(str(round(coefs[cursor])))
            cursor = cursor + 1
        elif current_token[0] == TOKEN_SPLIT:
            str_buf.append(" ---> ")
            if coefs[cursor] != 1:
                str_buf.append(str(round(coefs[cursor])))
            cursor = cursor + 1
    return "".join(str_buf)


def compute(input_str):
    lex = lexer(input_str)
    sorted_dict = extract_unique(lex)

    # print(lex)

    matrix = []
    prev_element = None
    positive = 1
    row = [0] * len(sorted_dict)
    for i, token in enumerate(lex):
        if token[0] == TOKEN_ELEMENT:
            prev_element = token[1]
            if i < len(lex) - 1 and lex[i + 1][0] != TOKEN_NUMBER:
                loc = sorted_dict.index(prev_element)
                row[loc] = positive
        elif token[0] == TOKEN_NUMBER:
            if prev_element is None:
                raise Exception("不是合法的化学反应式，发生了语法错误。")
            loc = sorted_dict.index(prev_element)
            row[loc] = positive * token[1]
        elif token[0] == TOKEN_PLUS:
            matrix.append(numpy.copy(row))
            row = [0] * len(sorted_dict)
            prev_element = None
        elif token[0] == TOKEN_SPLIT:
            matrix.append(numpy.copy(row))
            row = [0] * len(sorted_dict)
            positive = -1
            prev_element = None
    if lex[len(lex) - 1][0] != TOKEN_NUMBER:
        loc = sorted_dict.index(prev_element)
        row[loc] = positive
    matrix.append(numpy.copy(row))
    row = [0] * len(sorted_dict)

    if positive == 1:
        raise Exception("不是合法的化学反应式。")

    b = numpy.zeros(len(matrix))

    A = numpy.array(matrix).transpose()

    # print(numpy.array(matrix).transpose())
    # print("Rank A: {}, N: {}".format(matrix_rank(A), A.shape[0]))

    rank_A = matrix_rank(A)

    if rank_A == A.shape[1]:
        raise Exception(
            "该化学反应式无法配平，化学反应式输入不正确。齐次线性方程组有解的条件是系数矩阵的秩小于未知数的个数。可能是发生了语法或拼写错误，比如注意O(欧)和0(零)的区别。")

    solution = null_space(A)

    raw_x = solution[:, 0]

    x = numpy.copy(raw_x)

    # 求最大公约数，美化配平系数
    max_number = max(x)

    for i in range(len(x)):
        x[i] /= max_number

    # print(x)

    g = round(x[0] * 1000)

    for i in range(len(x)):
        g = math.gcd(g, round(x[i] * 1000))

    # print(g)

    for i in range(len(x)):
        x[i] *= (1 / g) * 1000

    # print(x)

    return format_result(lex, x)


if __name__ == "__main__":
    assert compute("C3H8 + O2 ---> CO2 + H2O") == "C3H8 + 5O2 ---> 3CO2 + 4H2O"
    assert compute("CH4 + O2 ---> CO2 + H2O") == "CH4 + 2O2 ---> CO2 + 2H2O"
    assert compute("NaHCO3 ---> Na2CO3 + H2O + CO2") == "2NaHCO3 ---> Na2CO3 + H2O + CO2"
    assert compute("CaCO3 ---> CaO + CO2") == "CaCO3 ---> CaO + CO2"
    assert compute("C + SiO2  --->  Si + CO") == "2C + SiO2 ---> Si + 2CO"

    print(compute(input()))

