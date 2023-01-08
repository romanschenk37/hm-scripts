import IPython.display as dp
import sympy as sp

from util import backwards_substitution, pivot


def gauss(a: sp.Matrix, b: sp.Matrix, pivoting: bool = False, output: bool = False) -> sp.Matrix:
    swaped_rows_counter = 0
    if output:
        dp.display(dp.Math(f'A = {sp.latex(a)}, \\quad b = {sp.latex(b)}'))
        dp.display(dp.Markdown('## Gauss'))

    n = len(b)
    u = a.copy()
    b = b.copy()

    for k in range(n - 1):
        if u[k, k].is_zero or pivoting:
            if(pivot(u, {'b': b}, k, output=output)):
                swaped_rows_counter += 1

        for i in range(k + 1, n):
            if not u[i, k].is_zero:
                factor = u[i, k] / u[k, k]
                u[i, k:n] -= factor * u[k, k:n]
                b[i] -= factor * b[k]

                if output:
                    dp.display(dp.Math(
                        f'z_{{{i + 1}}} \\equiv z_{{{i + 1}}} - ({sp.latex(factor)}) \\cdot z_{{{k + 1}}} \\Rightarrow ' + sp.latex(
                            u) + ' \\; ' + sp.latex(b)))

    if output:
        dp.display(dp.Markdown('## Determinante'))
        det_calc = ('det(A) = (-1)^' + str(swaped_rows_counter) + " * ")
        for i in range(n):
            det_calc += (' ' + str(u[i, i]) + ' ')
            if(i != n-1):
                det_calc += ' * '
            else:
                det_calc += (' = ' + str(a.det()))
        dp.display(dp.Math(det_calc))

        dp.display(dp.Markdown('## Rückwärtasdasdasfseinsetzen'))


    x = backwards_substitution(u, b, symbol='x', output=output)

    if output:
        dp.display(dp.Math('x = ' + sp.latex(x)))

    return x


# TODO fehlerfortpflanzung
