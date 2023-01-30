class WingParameters(object):
    def __init__(self):
        c_mean = 0.90
        self.b = c_mean/2
        self.tipo = 0 # 0: Wing; 1: Horizonal Tail; 2: Vertical Tail

