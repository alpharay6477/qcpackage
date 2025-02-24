
class MarcusCalculator:

    def __init__(self):
        self.HRindex = [2,1,8]
    def read_HRindex(self,HRindex):
        HRindex = [2,1,2]
        self.HRindex = HRindex
        return HRindex
    aon=321


mc = MarcusCalculator()
HRindex = []
hrindex2 = mc.read_HRindex(HRindex)
hc = mc.HRindex
print(hc)
print(hrindex2)
mc.HRindex[2] = 4

print(hrindex2)
print(mc.HRindex)