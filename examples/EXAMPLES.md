# Example potentials to try out

## Harmonic oscillator
ho = lambda x: 0.5 * x ** 2

## Square well
isw = lambda x, w, U: U * (abs(x) > w)

## Linear well
tp = lambda x, k: k * abs(x)

## Quartic oscillator
qo = lambda x: 0.5 * x ** 4

## Anharmonic oscillator
aho = lambda x: 0.25 * x ** 2 + 0.5 * x ** 4

## Double-well potential
dwp = lambda x, a: 0.5 * (x ** 2 - a ** 2) ** 2

## Linear ramp
lr = lambda x, k: k * x * (x > 0)

## Morse potential
mp = lambda x: 30 * (1 - np.exp(-1 * (x - 1))) ** 2

## Asymmetric double well potential
adwp = lambda x: 0.5 * ((x - 2) ** 2) * ((x + 2) ** 2 + 1)

# semi-circular potential
scp = lambda x, a: -np.sqrt(a ** 2 - x ** 2)