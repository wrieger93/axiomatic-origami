import itertools
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
        if type(self) != type(other):
            return False
        return self.x.equals(other.x) and self.y.equals(other.y)

    def __neq__(self, other):
        """Tests for inequality. Opposite of __eq__"""
        return not self.__eq__(other)

    def __hash__(self):
        """Returns a hash of the vector."""
        return hash((self.x, self.y))

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

    def __pos__(self):
        """Unary positive operator, doesn't change vector."""
        return self

    def __neg__(self):
        """Unary negation, equivalent to multiplying by -1."""
        return -1*self

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

    def project_onto_vector(self, other):
        """The vector projected onto the vector other."""
        return self.dot(other)/(other.dot(other)) * other

    def project_onto_line(self, line):
        """The point projected onto the line."""
        return (self-line.p).project_onto_vector(line.d) + line.p

    def lies_on_line(self, line):
        """True if the point lies on the line."""
        return self == self.project_onto_line(line)

    def lies_on_line_segment(self, lineseg):
        """True if the point ies on the line segment."""
        param = (self - lineseg.p1).dot(lineseg.p2 - lineseg.p1) / ((lineseg.p2 - lineseg.p1).dot(lineseg.p2 - lineseg.p1))
        return self.lies_on_line(lineseg.line_through()) and 0 <= param <= 1

    def line_dist(self, line):
        """The distance between the point and the line."""
        return (self - self.project_onto_line(line)).norm()

    def reflect_across(self, line):
        """The point reflected across the line."""
        return 2*self.project_onto_line(line) - self


class Line(object):
    """A class representing a two-dimensional line.
    
    The line is defined by a point on the line and a vector parallel to the line."""

    def __init__(self, p=None, d=None):
        """Initializes a line going through p parallel to d."""
        self.p = p
        self.d = d

        if self.p is None:
            self.p = Vector(0,0)
        if self.d is None:
            self.d = Vector(1,0)

    def __repr__(self):
        """Returns a string representation of the line."""
        return "({} + t*{})".format(self.p, self.d)

    def __eq__(self, other):
        """Tests if the two lines are equivalent representations of the same line."""
        if type(self) != type(other):
            return False
        return self.parallel_to(other) and self.p.lies_on_line(other)

    def __neq__(self, other):
        """Opposite of __eq__."""
        return not (self == other)

    def parallel_to(self, line):
        """True if the two lines are parallel."""
        return self.d.cross(line.d) == 0

    def reflect_across(self, line):
        """The given line reflected about the other line."""
        p1 = self.p.reflect_across(line)
        p2 = (self.p + self.d).reflect_across(line)
        return Line(p1, p2 - p1)

    def intersects_line(self, line):
        """The point where the two lines intersect. If the lines are parallel then
        None is returned. If the lines are the same then -1 is returned.
        """
        if self == line:
            return -1
        if self.parallel_to(line):
            return None

        A = sympy.Matrix([[self.d.x, -line.d.x],[self.d.y, -line.d.y]])
        b = sympy.Matrix([line.p.x - self.p.x, line.p.y - self.p.y])
        ans = A.LUsolve(b)
        return self.p + ans[0]*self.d

    def intersects_line_segment(self, lineseg):
        """The point where the line intersects the line segment. If the line
        segment lies on the line then -1 is returned. If the two don't intersect
        then None is returned.
        """
        if self == lineseg.line_through():
            return -1
        elif self.parallel_to(lineseg.line_through()):
            return None

        p = self.intersects_line(lineseg.line_through())
        if p.lies_on_line_segment(lineseg):
            return p
        else:
            return None


class LineSegment(object):
    """A class representing a two-dimensional line segment.

    The line segment is defined by the two endpoints of the segment."""

    def __init__(self, p1=None, p2=None):
        """Initializes a line segment with the endpoints p1 and p2."""
        self.p1 = p1
        self.p2 = p2

        if self.p1 is None:
            self.p1 = Vector(0,0)
        if self.p2 is None:
            self.p2 = Vector(1,0)

    def __repr__(self):
        """Returns a string representation of the line segment."""
        return "({} to {})".format(self.p1, self.p2)

    def __eq__(self, other):
        """Tests the line segments for equality."""
        if type(self) != type(other):
            return False
        return (self.p1 == other.p1 and self.p2 == other.p2) or (self.p1 == other.p2 and self.p2 == other.p1)

    def __neq__(self, other):
        """Opposite of __eq__."""
        return not (self == other)

    def __hash__(self):
        """Returns a hast of the line segment."""
        if self.p1.x < self.p2.x:
            return hash((self.p1, self.p2))
        else:
            return hash((self.p2, self.p1))

    def line_through(self):
        """The line passing through the line segment."""
        return Line(self.p1, self.p2 - self.p1)

    def intersects_line(self, line):
        """The point where the line segment intersects the line. Refer to the
        Line method intersects_line_segment()."""
        return line.intersects_line_segment(self)

    def intersects_line_segment(self, lineseg):
        """The point where the two line segments intersect. If they overlap at
        more than one point then -1 is returned. If they don't intersect then
        None is returned."""

        # if the line segments are parallel there's annoying cases
        if self.line_through().parallel_to(lineseg.line_through()):
            # if they're parallel and not on the same line they don't intersect
            if self.p1.lies_on_line(lineseg.line_through()):
                # literal edge cases
                # line segments are parallel but only endpoints coincide
                if self.p1 == lineseg.p1:
                    if (self.p2.lies_on_line_segment(lineseg) or
                            lineseg.p2.lies_on_line_segment(self)):
                        return -1
                    else:
                        return self.p1
                elif self.p1 == lineseg.p2:
                    if (self.p2.lies_on_line_segment(lineseg) or
                            lineseg.p1.lies_on_line_segment(self)):
                        return -1
                    else:
                        return self.p1
                elif self.p2 == lineseg.p1:
                    if (self.p1.lies_on_line_segment(lineseg) or
                            lineseg.p2.lies_on_line_segment(self)):
                        return -1
                    else:
                        return self.p2
                elif self.p2 == lineseg.p2:
                    if (self.p1.lies_on_line_segment(lineseg) or
                            lineseg.p1.lies_on_line_segment(self)):
                        return -1
                    else:
                        return self.p2
                else:
                    # none of the endpoints coincide so either they don't
                    # intersect or they overlap
                    if (self.p1.lies_on_line_segment(lineseg) or
                            self.p2.lies_on_line_segment(lineseg) or
                            lineseg.p1.lies_on_line_segment(self) or
                            lineseg.p2.lies_on_line_segment(sekf)):
                        return -1
                    else:
                        return None

            else:
                return None

        p = self.line_through().intersects_line(lineseg.line_through())
        if p.lies_on_line_segment(self) and p.lies_on_line_segment(lineseg):
            return p
        else:
            return None


class OrigamiPaper(object):
    """A class representing an arbitrary piece of origami paper.

    Axioms 1-7 can be applied to the paper with appropriate arguments.
    The class stores a list of points and line segments found. Only one line
    can be created per fold but the intersection of that line with every other 
    fold on the paper is added to the list.

    Possible future additions:
    * Allowing for an arbitrary paper shape (or maybe only convex)
    * Storing a history of folds made
    """
    
    def __init__(self):
        """Initialize the paper. This includes the boundary, the found points,
        and the found line segments."""
        self.points = [Vector(0,0), Vector(1,0), Vector(1,1), Vector(0,1)]
        self.linesegs = []
        for p in range(len(self.points)):
            self.linesegs.append(LineSegment(
                self.points[p],
                self.points[(p+1)%len(self.points)]))
        self.boundary = self.linesegs[:]

        self.points = set(self.points)
        self.linesegs = set(self.linesegs)

    def __repr__(self):
        """Returns a string representation of the paper."""
        return "boundary: {}\npoints: {}\nlinesegs: {}".format(
                self.boundary, self.points, self.linesegs)

    def intersects_boundary(self, line):
        """Returns a list of points where the given line intersects the
        boundary. Can be between 0 and 2 points. If the line coincides with
        one of the edges, -1 is returned."""
        points = []
        for seg in self.boundary:
            res = seg.intersects_line(line)
            if res == -1:
                return -1
            elif res is None:
                continue
            else:
                if res not in points:
                    points.append(res)

        return points

    def add_all_intersections(self, line):
        """Adds to self.points all the points where the given line
        intersects all other line segments. Also adds the line segment if
        applicable. If the line segment already exists
        or the line doesn't intersect the boundary in two points -1 is returned."""
        points = self.intersects_boundary(line)
        if points == -1 or len(points) < 2:
            return -1
        lineseg = LineSegment(points[0], points[1])

        if lineseg in self.linesegs:
            return -1

        for seg in self.linesegs:
            res = seg.intersects_line_segment(lineseg)
            if not (res == -1 or res is None):
                self.points.add(res)

        self.linesegs.add(lineseg)

    def axiom_1(self, p1, p2):
        """Returns the fold line through points p1 and p2. If p1 and p2 are
        for some reason equal None is returned. """
        if p1 == p2:
            return None
        return Line(p1, p2-p1)

    def axiom_2(self, p1, p2):
        """Returns the fold line that places p1 onto p2 (or vice versa).
        If p1 an p2 are for some reason equal None is returned."""
        if p1 == p2:
            return None
        midpoint = (p1 + p2)/sympy.S("2")
        return Line(midpoint, (p2-p1).rotate(sympy.pi/2))

    def axiom_4(self, p, seg):
        """Returns the fold line passing through point p perpendicular to l."""
        return Line(p, (seg.p2-seg.p1).rotate(sympy.pi/2))

if __name__ == "__main__":
    o = OrigamiPaper()
    for trial in range(3):
        tot = len(o.points)*(len(o.points)-1)/2
        count = 0
        combos = itertools.combinations(o.points, 2)
        for sub in combos:
            count += 1
            if count % 10 == 0:
                print(count, "out of", tot, sep=" ")
            p1, p2 = sub
            foldline = o.axiom_1(p1, p2) 
            if foldline is not None:
                o.add_all_intersections(foldline)
            foldline = o.axiom_2(p1, p2) 
            if foldline is not None:
                o.add_all_intersections(foldline)

