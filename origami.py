class Vector(object):
    def __init__(self, x=0, y=0):
        self.x = x
        self.y = y

    def __repr__(self):
        return "({}, {})".format(self.x, self.y)

    def __eq__(self, other):
        """Tests for equality."""
        return self.x == other.x and self.y == other.y

    def __neq__(self, other):
        """Tests for inequality. Opposite of __eq__"""
        return not self.__eq__(other)

    def __add__(self, other):
        """Adds two vectors."""
        return Vector(self.x + other.x, self.y + other.y)

    def __sub__(self, other):
        """Subtracts two vectors."""
        return Vector(self.x - other.x, self.y - other.y)

    def __mul__(self, other):
        """Multiplies a vector by a scalar."""
        return Vector(self.x*other, self.y*other)

    def __rmul__(self, other):
        """Multiplies a scalar by a vector."""
        return self * other 

    def __truediv__(self, other):
        """Divides a vector by a scalar."""
        return self*(1/other)

    def dot(self, other):
        """The dot product of two vectors."""
        return self.x * other.x + self.y * other.y

    def cross(self, other):
        """The magnitude of the 3D cross product between
        two vectors. Useful for testing if vectors are parallel.
        """
        return self.x * other.y - self.y * other.x

if __name__ == "__main__":
    a = Vector(1, 3)
    b = Vector(2, 5)
    c = 3
    print("a: {}, b: {}, c: {}".format(a, b, c))
    print("a+b:", a + b)
    print("a-b:", a - b)
    print("a*c:", a * c)
    print("c*a:", c * a)
    print("b/c:", b / c)
