from math import gcd

class FastFraction:

    def __init__(self, n, d):
        if d < 0:
            n = -n
            d = -d
        g = gcd(n, d)
        self.n = n // g
        self.d = d // g

    def __gt__(self, other):
        return self.n * other.d > other.n * self.d

    def __ge__(self, other):
        return self.n * other.d >= other.n * self.d

    def __lt__(self, other):
        return self.n * other.d < other.n * self.d

    def __le__(self, other):
        return self.n * other.d <= other.n * self.d

    def __eq__(self, other):
        return (self.n, self.d) == (other.n, other.d)

    def __neg__(self):
        return FastFraction(-self.n, self.d)

    def __mul__(self, other):
        return FastFraction(self.n * other.n, self.d * other.d)

    def __add__(self, other):
        return FastFraction(self.n * other.d + other.n * self.d, self.d * other.d)

    def __float__(self):
        return self.n / self.d

    def __str__(self):
        return str(float(self))

    def __repr__(self):
        return 'FastFraction({}, {})'.format(self.n, self.d)
