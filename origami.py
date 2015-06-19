from enum import Enum
import itertools
import logging
import sympy
import matplotlib.pyplot as plt

accuracy = 50

class IntersectionType(Enum):
    """An enum representing the possible intersections between two lines.

    Two lines (or line segments) in the X-Y plane can either
    1. intersect at exactly one point
    2. not intersect at all
    3. intersect at an infinite number of points

    This enum encapsulates these different possibilities to make dealing
    with line intersections easier.
    """

    # this can also be used to represent a finite number of intersection
    # points between a set of lines
    single = 1

    # no intersection point
    none = 2

    # lines overlap
    infinite = 3


class Vector(object):
    """A class representing a two-dimensional vector."""

    def __init__(self, x=0, y=0):
        """Initializes a point with x and y components.

        x and y should be either integers or a string representation of
        a number (e.g. "2/5" or "1/sqrt(2)") to make use of sympy's exact
        math capabilities.
        """
        # use sympy.S to use exact math
        self.x = sympy.S(x)
        self.y = sympy.S(y)

    def __repr__(self):
        """Returns a string representation of the vector."""
        return "<{}, {}>".format(self.x.round(6), self.y.round(6))

    def __eq__(self, other):
        """Tests for equality."""
        if type(self) != type(other):
            return False
        return (sympy.N(self.x, accuracy) == sympy.N(other.x, accuracy)
                and sympy.N(self.y, accuracy) == sympy.N(other.y, accuracy))

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
        if other.equals(0):
            raise ValueError("Cannot divide a vector by 0")
        return self*(1/other)

    def __pos__(self):
        """Unary positive operator, doesn't change vector."""
        return self

    def __neg__(self):
        """Unary negation, equivalent to multiplying by -1."""
        return -1*self

    def simplify(self):
        """Simplifies the vector's coordinates with sympy."""
        return Vector(sympy.simplify(self.x), sympy.simplify(self.y))

    def dot(self, other):
        """The dot product of two vectors."""
        return self.x * other.x + self.y * other.y

    def cross(self, other):
        """The magnitude of the 3D cross product between two vectors."""
        return self.x * other.y - self.y * other.x

    def norm(self):
        """The magnitude of the vector."""
        return sympy.sqrt(self.x * self.x + self.y * self.y)

    def normalize(self):
        """The unit vector in the same direction as the given vector."""
        return (self/self.norm()).simplify()

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

    def lies_on_line_segment(self, lseg):
        """True if the point lies on the line segment."""
        lseg_dir = lseg.p2 - lseg.p1
        param = (self - lseg.p1).dot(lseg_dir) / lseg_dir.dot(lseg_dir)
        return self.lies_on_line(lseg.line_through()) and 0 <= param <= 1

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
        """Tests if the two lines are nondistinct."""
        if type(self) != type(other):
            return False
        return self.parallel_to(other) and self.p.lies_on_line(other)

    def __neq__(self, other):
        """Opposite of __eq__."""
        return not (self == other)

    def simplify(self):
        """Simplifies the line using sympy."""
        return Line(self.p.simplify(), self.d.simplify())

    def parallel_to(self, line):
        """True if the two lines are parallel."""
        return sympy.N(sympy.Abs(self.d.cross(line.d)), accuracy) < sympy.S("10^-{}".format(accuracy))

    def reflect_across(self, line):
        """The given line reflected about the other line."""
        p1 = self.p.reflect_across(line)
        p2 = (self.p + self.d).reflect_across(line)
        return Line(p1, p2 - p1).simplify()

    def intersects_line(self, line):
        """The point where the two lines intersect.

        Returns the tuple (type, point), where type is an IntersectionType
        and point is the point (None if type is none or infinite).
        """
        if self == line:
            return IntersectionType.infinite, None
        if self.parallel_to(line):
            return IntersectionType.none, None

        A = sympy.Matrix([[self.d.x, -line.d.x],[self.d.y, -line.d.y]])
        b = sympy.Matrix([line.p.x - self.p.x, line.p.y - self.p.y])
        ans = A.LUsolve(b)
        return IntersectionType.single, self.p + ans[0]*self.d

    def intersects_line_segment(self, lineseg):
        """The point where the line intersects the line segment.

        Returns the tuple (type, point), where type is an IntersectionType
        and point is the point (None if type is none or infinite).
        """

        int_type, p = self.intersects_line(lineseg.line_through())
        if (int_type == IntersectionType.infinite or
                int_type == IntersectionType.none):
            return int_type, None
        elif p.lies_on_line_segment(lineseg):
            return IntersectionType.single, p
        else:
            return IntersectionType.none, None


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
        return ((self.p1 == other.p1 and self.p2 == other.p2) or
                (self.p1 == other.p2 and self.p2 == other.p1))

    def __neq__(self, other):
        """Opposite of __eq__."""
        return not (self == other)

    def __hash__(self):
        """Returns a hast of the line segment."""
        if self.p1.x < self.p2.x:
            return hash((self.p1, self.p2))
        else:
            return hash((self.p2, self.p1))

    def simplify(self):
        """Simplifies the line segment using sympy."""
        return LineSegment(self.p1.simplify(), self.p2.simplify())

    def line_through(self):
        """The line passing through the line segment."""
        return Line(self.p1, self.p2 - self.p1)

    def reflect_across(self, line):
        """The line segment reflected across the line."""
        return LineSegment(self.p1.reflect_across(line), self.p2.reflect_across(line)).simplify()

    def intersects_line(self, line):
        """The point where the line intersects the line segment.

        Just calls the corresponding method of the Line class for line
        segment intersection.

        Returns the tuple (type, point), where type is an IntersectionType
        and point is the point (None if type is none or infinite).
        """
        return line.intersects_line_segment(self)

    def intersects_line_segment(self, lseg):
        """The point where the two line segments intersect.

        Returns the tuple (type, point), where type is an IntersectionType
        and point is the point (None if type is none or infinite).
        """
        #  line segments are parallel
        if self.line_through().parallel_to(lseg.line_through()):
            # if they're parallel and on the same line there's some edge cases
            if self.line_through() == lseg.line_through():
                # line segments are parallel but only endpoints coincide
                if self.p1 == lseg.p1:
                    if (self.p2.lies_on_line_segment(lseg) or
                            lseg.p2.lies_on_line_segment(self)):
                        return IntersectionType.infinite, None
                    else:
                        return IntersectionType.single, self.p1
                elif self.p1 == lseg.p2:
                    if (self.p2.lies_on_line_segment(lseg) or
                            lseg.p1.lies_on_line_segment(self)):
                        return IntersectionType.infinite, None
                    else:
                        return IntersectionType.single, self.p1
                elif self.p2 == lseg.p1:
                    if (self.p1.lies_on_line_segment(lseg) or
                            lseg.p2.lies_on_line_segment(self)):
                        return IntersectionType.infinite, None
                    else:
                        return IntersectionType.single, self.p2
                elif self.p2 == lseg.p2:
                    if (self.p1.lies_on_line_segment(lseg) or
                            lseg.p1.lies_on_line_segment(self)):
                        return IntersectionType.infinite, None
                    else:
                        return IntersectionType.single, self.p2
                else:
                    # none of the endpoints coincide so either they don't
                    # intersect or they overlap
                    if (self.p1.lies_on_line_segment(lseg) or
                            self.p2.lies_on_line_segment(lseg) or
                            lseg.p1.lies_on_line_segment(self) or
                            lseg.p2.lies_on_line_segment(sekf)):
                        return IntersectionType.infinite, None
                    else:
                        return IntersectionType.none, None

            else:
                # if the lines are parallel and aren't the same,
                # they don't intersect
                return IntersectionType.none, None

        # line segments aren't parallel
        else:
            int_type, p = self.line_through().intersects_line(lseg.line_through())
            # in theory this should always be false, but here for safety
            if int_type != IntersectionType.single:
                return int_type, None
            elif p.lies_on_line_segment(self) and p.lies_on_line_segment(lseg):
                return IntersectionType.single, p
            else:
                return IntersectionType.none, None


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

        self.points = self.points
        self.linesegs = self.linesegs

    def __repr__(self):
        """Returns a string representation of the paper."""
        return "boundary: {}\npoints: {}\nlinesegs: {}".format(
                self.boundary, self.points, self.linesegs)

    def draw(self):
        for lseg in self.boundary:
            plt.plot([lseg.p1.x, lseg.p2.x], [lseg.p1.y, lseg.p2.y], "k-")
        for lseg in self.linesegs:
            plt.plot([lseg.p1.x, lseg.p2.x], [lseg.p1.y, lseg.p2.y], "k-")
        plt.show()

    def intersects_boundary(self, line):
        """A list of points where the given line intersects the boundary.

        Returns a tuple (intersection_type, points), where points is a list
        of intersection points. If intersection_type is infinite, the line
        overlaps one of the boundary lines. If it's none, the line doesn't
        intersect the boundary. Otherwise, points will be a list of length
        one or two.
        """
        points = []
        for seg in self.boundary:
            intersection_type, point = seg.intersects_line(line)
            if intersection_type == IntersectionType.infinite:
                return IntersectionType.infinite, None
            elif intersection_type == IntersectionType.none:
                continue
            else:
                if point not in points:
                    points.append(point)

        return IntersectionType.single, points

    def add_all_intersections(self, line):
        """Adds to self.points all the points where the given line
        intersects all other line segments. Also adds the line segment if
        applicable.
        """
        intersection_type, points = self.intersects_boundary(line)
        if intersection_type != IntersectionType.single or len(points) < 2:
            return
        lineseg = LineSegment(points[0], points[1])

        if lineseg in self.linesegs:
            return

        self.linesegs.append(lineseg)

        for seg in self.linesegs:
            intersection_type, point = seg.intersects_line_segment(lineseg)
            if intersection_type == IntersectionType.single and point not in self.points:
                self.points.append(point)

    def axiom_1(self, p1, p2):
        """Returns the fold line through points p1 and p2."""
        if p1 == p2:
            raise ValueError("Points are not distinct.")
        return Line(p1, p2-p1)

    def axiom_2(self, p1, p2):
        """Returns the fold line that places p1 onto p2 (or vice versa)."""
        if p1 == p2:
            raise ValueError("Points are not distinct.")
        midpoint = (p1 + p2)/sympy.S("2")
        return Line(midpoint, (p2-p1).rotate(sympy.pi/2))

    def axiom_3(self, lseg1, lseg2):
        """Returns a list of fold lines that place lseg1 onto lseg2.

        There can be up to two folds that accomplish this. The function makes
        sure that the reflected segments overlap (or else you couldn't
        align the fold).
        """
        if lseg1 == lseg2:
            raise ValueError("Line segments are not distinct.")
        int_type, p = lseg1.intersects_line_segment(lseg2)
        if int_type == IntersectionType.infinite:
            return []
        elif lseg1.line_through().parallel_to(lseg2.line_through()):
            fold_line = Line((lseg1.p1 + lseg2.p1)/sympy.S(2), lseg1.line_through().d)
            int_type, p =  lseg1.reflect_across(fold_line).intersects_line_segment(lseg2)
            if int_type != IntersectionType.infinite:
                return []
            else:
                return [fold_line]
        else:
            _, int_point = lseg1.line_through().intersects_line(lseg2.line_through())
            u1 = lseg1.line_through().d.normalize()
            u2 = lseg2.line_through().d.normalize()
            lst = []
            fold1 = Line(int_point, u1 + u2)
            fold2 = Line(int_point, u1 - u2)
            int_type1, p1 =  lseg1.reflect_across(fold1).intersects_line_segment(lseg2)
            int_type2, p2 =  lseg1.reflect_across(fold2).intersects_line_segment(lseg2)
            if int_type1 == IntersectionType.infinite:
                lst.append(fold1)
            if int_type2 == IntersectionType.infinite:
                lst.append(fold2)
            return lst

    def axiom_4(self, p, seg):
        """Returns the fold line passing through point p perpendicular to l."""
        return Line(p, (seg.p2-seg.p1).rotate(sympy.pi/2))


class TerminalFolder(object):
    """A class to allow interactive manipulation of the OrigamiPaper class."""

    def __init__(self, paper):
        """Initialize the class with an OrigamiPaper instance."""
        self.paper = paper

    def get_input(self, choices_dict, prompt_msg, error_msg="Try again.", print_choices=False):
        """Prompts the user for input from a set of choices.

        choices_dict should have entries of the form (string: obj) since
        input() returns a string.
        """
        if print_choices:
            print("Choices:")
            for k, v in sorted(choices_dict.items()):
                print("{}: {}".format(k, v))
            print()
        res = None
        while True:
            res = input(prompt_msg + " ")
            if res not in choices_dict:
                print(error_msg)
            else:
                break
        print()
        return choices_dict[res]

    def fold_interactive(self):
        """An interactive prompt to test folding.

        Allows for terminal input to try out the different folding axioms.
        """
        axioms_dict = dict([(str(i), i) for i in range(1, 8)])
        done = False
        while not done:
            self.paper.draw()

            points_dict = dict([(str(i), p) for i, p in enumerate(self.paper.points)])
            linesegs_dict = dict([(str(i), lseg) for i, lseg in enumerate(self.paper.linesegs)])

            axiom_num = self.get_input(axioms_dict, "Axiom? [1-7]")

            if axiom_num == 1:
                p1 = self.get_input(points_dict, "First point?", print_choices=True)
                p2 = self.get_input(points_dict, "Second point?", print_choices=False)
                self.paper.add_all_intersections(self.paper.axiom_1(p1, p2))

            elif axiom_num == 2:
                p1 = self.get_input(points_dict, "First point?", print_choices=True)
                p2 = self.get_input(points_dict, "Second point?", print_choices=False)
                self.paper.add_all_intersections(self.paper.axiom_2(p1, p2))

            elif axiom_num == 3:
                lseg1 = self.get_input(linesegs_dict, "First line segment?", print_choices=True)
                lseg2 = self.get_input(linesegs_dict, "Second line segment?", print_choices=False)
                fold_lst = self.paper.axiom_3(lseg1, lseg2)
                if len(fold_lst) == 0:
                    print("No fold found.")
                elif len(fold_lst) == 1:
                    self.paper.add_all_intersections(fold_lst[0])
                else:
                    folds_dict = dict([(str(i), f) for i, f in enumerate(fold_lst)])
                    fold = self.get_input(folds_dict, "Which fold?", print_choices=True)
                    self.paper.add_all_intersections(fold)

            elif axiom_num == 4:
                p = self.get_input(points_dict, "Point?", print_choices=True)
                lseg = self.get_input(linesegs_dict, "Line segment?", print_choices=True)
                self.paper.add_all_intersections(self.paper.axiom_4(p, lseg))

            done = not self.get_input({"y": True, "Y": True, "n": False, "N": False}, "Continue? [y/n]")

        print("Goodbye!")


if __name__ == "__main__":
    o = OrigamiPaper()
    t = TerminalFolder(o)
    t.fold_interactive()
