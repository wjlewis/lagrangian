// type Num = number | [Num, Num];

export function plus(x, y) {
  if (typeof x === 'number') {
    if (typeof y === 'number') {
      return x + y;
    } else {
      return [plus(x, y[0]), y[1]];
    }
  } else {
    if (typeof y === 'number') {
      return [plus(x[0], y), x[1]];
    } else {
      const [a, b] = x;
      const [c, d] = y;
      return [plus(a, c), plus(b, d)];
    }
  }
}

export function times(x, y) {
  if (typeof x === 'number') {
    if (typeof y === 'number') {
      return x * y;
    } else {
      return [times(x, y[0]), times(x, y[1])];
    }
  } else {
    if (typeof y === 'number') {
      return [times(x[0], y), times(x[1], y)];
    } else {
      const [a, b] = x;
      const [c, d] = y;
      return [times(a, c), plus(times(b, c), times(a, d))];
    }
  }
}

export function neg(x) {
  if (typeof x === 'number') {
    return -x;
  } else {
    return [neg(x[0]), neg(x[1])];
  }
}

export function minus(x, y) {
  return plus(x, neg(y));
}

export function inv(x) {
  if (typeof x === 'number') {
    return 1 / x;
  } else {
    const [a, b] = x;
    return [inv(a), times(neg(b), inv(pow(a, 2)))];
  }
}

export function div(x, y) {
  return times(x, inv(y));
}

export function pow(x, n) {
  if (typeof x === 'number') {
    return x ** n;
  } else {
    const [a, b] = x;
    // (a + be)^n = a^n + na^{n-1}be
    return [pow(a, n), times(times(n, pow(a, n - 1)), b)];
  }
}

export function sin(x) {
  if (typeof x === 'number') {
    return Math.sin(x);
  } else {
    const [a, b] = x;
    // f(a + be) = f(a) + Df(a)be
    return [sin(a), times(cos(a), b)];
  }
}

export function cos(x) {
  if (typeof x === 'number') {
    return Math.cos(x);
  } else {
    const [a, b] = x;
    return [cos(a), times(neg(sin(a)), b)];
  }
}

export function exp(x) {
  if (typeof x === 'number') {
    return Math.exp(x);
  } else {
    const [a, b] = x;
    return [exp(a), times(exp(a), b)];
  }
}

export function eq(x, y) {
  if (typeof x === 'number') {
    if (typeof y === 'number') {
      return x === y;
    } else {
      return eq(q, y[0]);
    }
  } else {
    if (typeof y === 'number') {
      return eq(x[0], y);
    } else {
      return eq(x[0], y[0]);
    }
  }
}

export function lt(x, y) {
  if (typeof x === 'number') {
    if (typeof y === 'number') {
      return x < y;
    } else {
      return lt(q, y[0]);
    }
  } else {
    if (typeof y === 'number') {
      return lt(x[0], y);
    } else {
      return lt(x[0], y[0]);
    }
  }
}

export function lte(x, y) {
  return lt(x, y) || eq(x, y);
}

export function gt(x, y) {
  return !lte(x, y);
}

export function gte(x, y) {
  return gt(x, y) || eq(x, y);
}

// Differentiation

export function D(f) {
  return x => {
    const y = f([x, 1]);
    return y[1];
  };
}

// Lagrange

// xs : Num[]
// ys : Num[]
export function lagrangeInterpolation(xs, ys) {
  return x => {
    return ys.reduce((sum, y, i) => {
      const li = lagrangeBasis(xs, i);

      return plus(sum, times(y, li(x)));
    }, 0);
  };
}

function lagrangeBasis(xs, i) {
  const xi = xs[i];

  return x => {
    return xs.reduce((prod, xj, j) => {
      if (i === j) {
        return prod;
      } else {
        return times(prod, div(minus(x, xj), minus(xi, xj)));
      }
    }, 1);
  };
}

export function flatLagrangeInterpolation(xs, ys) {
  return x => {
    return ys.reduce((sum, yi, i) => {
      const li = flatLagrangeBasis(xs, i);

      return plus(sum, times(yi, li(x)));
    }, 0);
  };
}

export function flatLagrangeBasis(xs, i) {
  const xi = xs[i];
  const xStar = getXStar(xs, i);
  const li = lagrangeBasis(xs, i);

  return x => {
    return times(li(x), div(minus(x, xStar), minus(xi, xStar)));
  };
}

function getXStar(xs, i) {
  if (i === 0) {
    const x0 = xs[0];
    const n = prod(xs.slice(1).map(xj => minus(x0, xj)));

    const d = xs.reduce((sum, _xj, j) => {
      if (j === 0) {
        return sum;
      } else {
        const prod = xs.reduce((prod, xk, k) => {
          if (k === 0 || k === j) {
            return prod;
          } else {
            return times(prod, minus(x0, xk));
          }
        }, 1);

        return plus(sum, prod);
      }
    }, 0);

    return plus(x0, div(n, d));
  } else {
    return xs[0];
  }
}

export function sum(xs) {
  return xs.reduce(plus, 0);
}

export function prod(xs) {
  return xs.reduce(times, 1);
}

// Simpson's Rule

// x0, x1 : number
export function simpsonIntegral(f, x0, x1, step = 0.1) {
  if (x1 < x0) {
    return -simpsonIntegral(f, x1, x0, step);
  }

  const nSteps = Math.ceil((x1 - x0) / step);
  const st = (x1 - x0) / nSteps;

  let sum = 0;
  for (let i = 0; i < nSteps; i++) {
    const a = x0 + i * st;
    const b = a + st;
    const h = (a + b) / 2;

    sum += ((b - a) / 6) * (f(a) + 4 * f(h) + f(b));
  }

  return sum;
}

// World

export function compose(f, g) {
  return x => f(g(x));
}

// Transform a coordinate path into a function from time to the local tuple.
export function gamma(q) {
  const Dq = D(q);
  return t => [t, q(t), Dq(t)];
}
