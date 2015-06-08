import sympy

class Vector(object):
    """A class representing a two-dimensional vector."""
    
    def __init__(self, x=0, y=0):
        """Initializes a point with x and y components."""
        self.x = sympy.S(x)
        self.y = sympy.S(y)

    def __repr__(self):
        """Returns a string representation of the vector."""
        return "<{}, {}>".format(self.x, self.y)

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
        if other == 0:
            raise ValueError("Cannot divide a vector by 0")
        return self*(1/other)

    def dot(self, other):
        """The dot product of two vectors."""
        return self.x * other.x + self.y * other.y

    def cross(self, other):
        """The magnitude of the 3D cross product between
        two vectors. Useful for testing if vectors are parallel.
        """
        return self.x * other.y - self.y * other.x

    def norm(self):
        """The magnitude of the vector."""
        return sympy.sqrt(self.x * self.x + self.y * self.y)

    def normalize(self):
        """The unit vector in the same direction as the given vector."""
        return self/self.norm()

    def rotate(self, angle):
        """The vector rotated by angle radians counterclockwise."""
        return Vector(sympy.cos(angle)*self.x - sympy.sin(angle)*self.y,
                sympy.sin(angle)*self.x + sympy.cos(angle)*self.y)

    def project_onto(self, other):
        """The vector projected onto the vector other."""
        return self.dot(other)/(other.dot(other)) * other

if __name__ == "__main__":
    # check that all the operations work

    a = Vector(sympy.S(1), sympy.S(3))
    b = Vector(sympy.S(2), sympy.S(5))
    c = sympy.S(3)

    print("a: {}, b: {}, c: {}".format(a, b, c))
    print("a+b:", a + b)
    print("a-b:", a - b)
    print("a*c:", a * c)
    print("c*a:", c * a)
    print("b/c:", b / c)
    print("a.b:", a.dot(b))
    print("axb:", a.cross(b))
    print("bxa:", b.cross(a))
    print("|a|:", a.norm())
    print("normalized a:", a.normalize())
    print("a rotated by pi/2:", a.rotate(sympy.pi/2))
    print("a projected onto b:", a.project_onto(b))
    print("b projected onto a:", b.project_onto(a))
