from origami import *


class TerminalFolder(object):
    """A class to allow interactive manipulation of the OrigamiPaper class.

    From the name, the manipulation is done via a command prompt."""

    def __init__(self, paper):
        """Initialize the class with an OrigamiPaper instance."""
        self.paper = paper

    def draw_paper(self):
        plt.ion()
        plt.cla()
        plt.axis([0,1,0,1])
        # for lseg in self.boundary:
            # plt.plot([lseg.p1.x, lseg.p2.x], [lseg.p1.y, lseg.p2.y], "k-")
        for lseg in self.paper.linesegs:
            plt.plot([lseg.p1.x, lseg.p2.x], [lseg.p1.y, lseg.p2.y], "k-")
        ax = plt.gca()
        ax.set_aspect("equal")
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        plt.show(block=False)

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
            self.draw_paper()

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
