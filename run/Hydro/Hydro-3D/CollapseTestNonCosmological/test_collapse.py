from yt.mods import *
from yt.testing import *
from yt.utilities.answer_testing.api import AnswerTestingTest
from yt.utilities.answer_testing.framework import \
    requires_outputlog, \
    sim_dir_load

_pf_name = "CollapseTestNonCosmological.enzo"
_dir_name = os.path.dirname(__file__)

class TestCollapseMaxValue(AnswerTestingTest):
    _type_name = "MaxValue"
    _attrs = ()

    def __init__(self, sim):
        self.pf = sim
    
    def run(self):
        result = []
        for my_pf in self.pf:
            result.append(my_pf.h.find_max("Density")[0])
        return result

    def compare(self, new_result, old_result):
        for i in range(len(new_result)):
            assert_rel_equal(new_result[i], old_result[i], 2)

@requires_outputlog(_dir_name, _pf_name)
def test_collapse_max_value():
    sim = sim_dir_load(_pf_name, path=_dir_name, 
                       find_outputs=True)
    sim.get_time_series()
    
    yield TestCollapseMaxValue(sim)
