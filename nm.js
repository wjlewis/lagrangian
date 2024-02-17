/**
 * Nelder-Mead Downhill Simplex Method
 */
export function nelderMead(
  f,
  simplex,
  {
    alpha = 1,
    gamma = 2,
    rho = 1 / 2,
    sigma = 1 / 2,
    epsilon = 1e-12,
    maxIter = 10_000,
  } = {}
) {
  if (simplex.length < 3) {
    throw new Error('dimension must be at least 2');
  }

  const vs = simplex.map(x => ({ x, z: f(...x) }));

  for (let i = 0; i < maxIter; i++) {
    // 1. Order points.
    vs.sort((v1, v2) => v1.z - v2.z);
    if (stdDev(...vs.map(({ z }) => z)) < epsilon) {
      break;
    }

    // 2. Calculate the centroid of all but the worst point.
    const x0 = centroid(...vs.slice(0, -1).map(({ x }) => x));

    // Worst, second-worst, and best points and values.
    const { x: xW, z: zW } = vs[vs.length - 1];
    const { z: zS } = vs[vs.length - 2];
    const { x: xB, z: zB } = vs[0];

    // 3. Reflection.
    const xR = plus(x0, scale(alpha, minus(x0, xW)));
    const zR = f(...xR);

    if (zR < zS) {
      // Reflected point is better than the second worst.
      if (zR >= zB) {
        // Reflected point is not better than the best.
        vs[vs.length - 1] = { x: xR, z: zR };
        continue;
      } else {
        // Reflected point is the best so far.
        // 4. Expansion.
        const xE = plus(x0, scale(gamma, minus(xR, x0)));
        const zE = f(...xE);

        if (zE < zR) {
          // Expanded point is better than reflected point.
          vs[vs.length - 1] = { x: xE, z: zE };
          continue;
        } else {
          vs[vs.length - 1] = { x: xR, z: zR };
          continue;
        }
      }
    } else {
      // Reflected point is not better than the second worst.
      // 5. Contraction.
      if (zR < zW) {
        // Reflected point is better than worst point.
        const xC = plus(x0, scale(rho, minus(xR, x0)));
        const zC = f(...xC);

        if (zC < zR) {
          vs[vs.length - 1] = { x: xC, z: zC };
          continue;
        }
      } else {
        // Reflected point is the worst so far.
        const xC = plus(x0, scale(rho, minus(xW, x0)));
        const zC = f(...xC);

        if (zC < zW) {
          vs[vs.length - 1] = { x: xC, z: zC };
          continue;
        }
      }
    }

    // 6. Shrinking.
    for (let i = 1; i < vs.length; i++) {
      const xI = plus(xB, scale(sigma, minus(vs[i].x, xB)));
      const zI = f(...xI);
      vs[i] = { x: xI, z: zI };
    }
  }

  // Return the best point.
  return vs[0].x;
}

function stdDev(...xs) {
  return Math.sqrt(variance(...xs));
}

function variance(...xs) {
  const mu = mean(...xs);
  return mean(...xs.map(x => (x - mu) ** 2));
}

function mean(...xs) {
  const sum = xs.reduce((sum, x) => sum + x, 0);
  return sum / xs.length;
}

// Vector operations.

function centroid(...xs) {
  return scale(1 / xs.length, sum(...xs));
}

function scale(f, x) {
  if (typeof x === 'number') {
    return f * x;
  } else {
    return x.map(xi => scale(f, xi));
  }
}

function sum(...xs) {
  if (xs.length === 0) {
    throw new Error('cannot sum 0 points');
  }

  const [x, ...rest] = xs;
  return rest.reduce(plus, x);
}

function plus(x1, x2) {
  if (typeof x1 === 'number' && typeof x2 === 'number') {
    return x1 + x2;
  } else {
    return zipWith(x1, x2, (xi, yi) => plus(xi, yi));
  }
}

function minus(x1, x2) {
  return plus(x1, scale(-1, x2));
}

function zipWith(xs, ys, f) {
  const out = new Array(xs.length);
  for (let i = 0; i < xs.length; i++) {
    out[i] = f(xs[i], ys[i]);
  }
  return out;
}

export function lerp(x1, x2, n) {
  const out = new Array(n + 2);
  out[0] = x1;
  out[out.length - 1] = x2;
  const step = scale(1 / (out.length - 1), minus(x2, x1));

  for (let i = 1; i < out.length - 1; i++) {
    out[i] = plus(out[i - 1], step);
  }

  return out;
}
