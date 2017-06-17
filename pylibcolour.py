# See LICENSE file for copyright and license details.

import math as _math

def _pow(x, p):
    return -((-x) ** p) if x < 0 else x ** p

def _cbrt(x):
    return _pow(x, 1. / 3.)

def _transpose(M):
    return [[r[c] for r in M] for c in range(len(M[0]))]

def _multiply(A, B):
    return [[sum(x * y for x, y in zip(a, b)) for b in _transpose(B)] for a in A]

class Colour(object):
    __matrices = {} ## TODO

    def __init__(self):
        pass

    def delta_e(self, other):
        x, y = CIELAB(self), CIELAB(other)
        L = x.L - y.L
        a = x.a - y.a
        b = x.b - y.b
        L *= L
        a *= a
        b *= b
        return (L + a + b) ** 0.5

    def cache(self, other):
        try:
            if isinstance(self, RGB) or isinstance(other, RGB):
                to_ciexyz = self.M     if isinstance(self,  RGB) else Colour.__get_matrix(self.name, 'CIEXYZ')
                fr_ciexyz = other.Minv if isinstance(other, RGB) else Colour.__get_matrix('CIEXYZ', other.name)
                Colour.__matrices[(self.matrix_id(), other.matrix_id())] = _multiply(fr_ciexyz, to_ciexyz)
                fr_ciexyz = self.Minv if isinstance(self,  RGB) else Colour.__get_matrix('CIEXYZ', self.name)
                to_ciexyz = other.M   if isinstance(other, RGB) else Colour.__get_matrix(other.name, 'CIEXYZ')
                Colour.__matrices[(other.matrix_id(), self.matrix_id())] = _multiply(fr_ciexyz, to_ciexyz)
        except:
            pass

    def decache(self, other):
        if isinstance(self, RGB) or isinstance(other, RGB):
            try:
                del Colour.__matrices[(self.matrix_id(), other.matrix_id())]
            except:
                pass
            try:
                del Colour.__matrices[(other.matrix_id(), self.matrix_id())]
            except:
                pass

    def get_configuration(self):
        return {}

    def copy(self):
        return type(self)(*(self.get_params()), **(self.get_configuration()))

    def get_params(self, linear = None):
        if linear is None or linear:
            return self.get_linear_params()
        else:
            raise Exception('linear must be None or True')

    def set_params(self, *args, linear = None):
        if linear is None or linear:
            return self.set_linear_params(*args)
        else:
            raise Exception('linear must be None or True')

    def matrix_id(self):
        return self.name

    def __init(self, args, name, linear):
        self.name = name
        self.linear = linear
        if len(args) == 0:
            return (None, None, None)
        if len(args) == 1:
            args = args[0]
            if isinstance(args, Colour):
                return Colour.__convert(args, self)
        def conv(a):
            try:
                return float(a):
            except ValueError:
                return float.fromhex(a)
        return tuple(conv(a) for a in args)

    def __repr__(self):
        def r(x):
            if type(a) is float:
                return repr(a.hex())
            elif type(a) is tuple:
                if len(a) == 1:
                    return '(%r,)' % a[0]
                return '(%s)' % ', '.join(map(r, a))
            elif type(a) is list:
                return '[%s]' % ', '.join(map(r, a))
            else:
                return repr(a)
        va = self.get_params()
        kw = self.get_configuration()
        s = []
        for v in va:
            s.append(r(v))
        for k, v in kw:
            s.append('%s = %s' % (k, r(v)))
        return '%s(%s)' % (self.name, ', '.join(s))

    @staticmethod
    def __get_matrix(fr, to):
        return Colour.__matrices[(fr.matrix_id(), to.matrix_id())]

    @staticmethod
    def __matrix_convert(fr, to, a, b, c):
        M = Colour.__get_matrix(fr, to)
        return (M[1][1] * a + M[1][2] * b + M[1][3] * c,
                M[2][1] * a + M[2][2] * b + M[2][3] * c,
                M[3][1] * a + M[3][2] * b + M[3][3] * c)

    @staticmethod
    def __convert(fr, to):
        try:
            M = Colour.__get_matrix(fr, to)
            if M is not None:
                (a, b, c) = fr.get_params(linear = True)
                (a, b, c) = (M[1][1] * a + M[1][2] * b + M[1][3] * c,
                             M[2][1] * a + M[2][2] * b + M[2][3] * c,
                             M[3][1] * a + M[3][2] * b + M[3][3] * c)
                to = to.copy()
                to.set_params(a, b, c, linear = True)
                return to.get_params()
        except:
            pass

        if isinstance(to, RGB):
            have_transfer = False
            with_transfer = to.with_transfer
            while True:
                if isinstance(fr, RGB) and fr.M == to.M:
                    (R, G, B) = fr.get_params()
                    have_transfer = fr.with_transfer
                    if have_transfer and with_transfer and not to.transfer_function.same(fr.transfer_function):
                        (R, G, B) = fr.transfer_function.decode(R, G, B)
                        have_transfer = False
                    break
                if not isinstance(fr, CIEXYZ):
                    fr = CIEXYZ(fr)
                R = to.Minv[0][0] * fr.X + to.Minv[0][1] * fr.Y + to.Minv[0][2] * fr.Z
                G = to.Minv[1][0] * fr.X + to.Minv[1][1] * fr.Y + to.Minv[1][2] * fr.Z
                B = to.Minv[2][0] * fr.X + to.Minv[2][1] * fr.Y + to.Minv[2][2] * fr.Z
            if have_transfer != with_transfer:
                if with_transfer:
                    if to.transfer_function is not None:
                        (R, G, B) = to.transfer_function.encode(R, G, B)
                else:
                    if fr.transfer_function is not None:
                        (R, G, B) = fr.transfer_function.decode(R, G, B)
            return (R, G, B)

        elif isinstance(to, sRGB) and isinstance(fr, sRGB):
            if to.with_transfer == fr.with_transfer:
                return fr.get_params()
            else:
                return fr.get_params(linear = not to.with_transfer)

        elif isinstance(to, CIEUVW):
            if isinstance(fr, CIEUVW):
                w = fr.W * 13
                return (fr.U + w * (fr.u0 - to.u0), fr.V + w * (fr.v0 - to.v0), fr.W)
            else:
                if not isinstance(fr, CIE1960UCS):
                    fr = CIE1960UCS(fr)
                Y = 25. * _cbrt(fr.Y) - 17.
                w = Y * 13.
                return (w * (fr.u - to.u0), w * (fr.v - to.v0), fr.Y)

        elif isinstance(to, CIELUV):
            if isinstance(fr, CIELUhuv):
                pi2 = 2. * _math.pi
                if to.white.X != fr.white.X or to.white.Y != fr.white.Y or to.white.Z != fr.white.Z:
                    fr = CIELUhuv(CIEXYZ(fr), one_revolution = pi2)
                elif fr.one_revolution != pi2:
                    fr = CIELUhuv(fr.L, fr.C, fr.h * pi2 / fr.one_revolution, one_revolution = pi2)
                return (fr.L, fr.C * _math.cos(fr.h), fr.C * _math.sin(fr.h))
            else:
                if isinstance(fr, CIELUV):
                    if fr.white.X == to.white.X and fr.white.Y == to.white.Y and fr.white.Z == to.white.Z:
                        return fr.get_params()
                if not isinstance(fr, CIEXYZ):
                    fr = CIEXYZ(fr)
                wx, wy, x, y = to.white.X, to.white.Y, fr.x, fr.y
                wt = wx + 15. * wy + 3. * to.white.Z
                t  =  x + 15. *  y + 3. * fr.z
                u = 4. * (x / t - wx / wt)
                v = 9. * (y / t - wy / wt)
                L = y / wy
                L2 = L * 24389.
                L = L2 / 27. if L2 <= 216. else _cbrt(L) * 116. - 16.
                L2 = L * 13.
                return (L, u * L2, v * L2)

        elif isinstance(to, CIELChuv):
            if isinstance(fr, CIELChuv):
                if to.white.X == fr.white.X and to.white.Y == fr.white.Y and to.white.Z == fr.white.Z:
                    if to.one_revolution != fr.one_revolution:
                        h = fr.h * to.one_revolution / fr.one_revolution
                        fr = CIELChuv(fr.L, fr.C, h, white = fr.white, one_revolution = to.one_revolution)
                    return fr.get_params()
            if not isinstance(fr, CIELUV):
                fr = CIELUV(fr, white = to.white)
            elif to.white.X != fr.white.X or to.white.Y != fr.white.Y or to.white.Z != fr.white.Z:
                fr = CIELUV(fr, white = to.white)
            h = _math.atan2(fr.v, fr.u) / (2. * _math.pi) * to.one_revolution
            if h < 0:
                h += to.one_revolution
            return (fr.L, (fr.u * fr.u + fr.v * fr.v) ** 0.5, h)

        elif isinstance(fr, RGB) and isinstance(to, CIEXYZ):
            (R, G, B) = fr.get_params(linear = True)
            X = fr.M[0][0] * R + fr.M[0][1] * G + fr.M[0][2] * B
            Y = fr.M[1][0] * R + fr.M[1][1] * G + fr.M[1][2] * B
            Z = fr.M[2][0] * R + fr.M[2][1] * G + fr.M[2][2] * B
            return (X, Y, Z)

        elif to.name == fr.name:
            return fr.get_params()

        elif to.linear:
            if isinstance(fr, sRGB):
                if fr.with_transfer:
                    fr = sRGB(fr, with_transfer = False)
            elif not fr.linear:
                if not isinstance(to, CIEXYZ):
                    fr = CIEXYZ(fr)
                elif isinstance(fr, CIExyY):
                    if fr.y == 0:
                        return (fr.Y, fr.Y, fr.Y)
                    else:
                        return (fr.x * fr.Y / fr.y, fr.Y, (1. - fr.x - fr.y) * fr.Y / fr.y)
                elif isinstance(fr, CIELAB):
                    Y = CIELAB.decode_transfer((fr.L + 16.) / 116.)
                    X = CIELAB.decode_transfer((Y + fr.a) / 500.) * 0.95047
                    Z = CIELAB.decode_transfer((Y - fr.b) / 200.) * 1.08883
                    return (X, Y, Z)
                elif isinstance(fr, CIELChuv) or isinstance(fr, CIELUV):
                    if isinstance(fr, CIELChuv):
                        fr = CIELUV(fr, white = fr.white)
                    X = fr.white.X
                    Y = fr.white.Y
                    L13 = fr.L * 13.
                    t = X + 15. * Y + 3. * fr.white.Z
                    u = fr.U / L13 + 4 * X / t
                    v = fr.V / L13 + 9 * Y / t
                    if fr.L <= 8.:
                        Y *= fr.L * 27. / 24389.
                    else:
                        L = (fr.L + 16.) / 116.
                        Y *= L * L * L
                    X = 2.25 * Y * u / v
                    Z = Y * (3. / v - 0.75 * u / v - 5.)
                    return (X, Y, Z)
                elif isinstance(fr, CIEUVW) or isinstance(fr, CIE1960UCS):
                    if isinstance(fr, CIEUVW):
                        fr = CIE1960UCS(fr)
                    Y = fr.Y
                    X = 1.5 * Y * fr.u / fr.v
                    Z = (4. * Y - Y / fr.u - 10. * Y * fr.v) / (2 * fr.v)
                    return (X, Y, Z)
                else:
                    fr = sRGB(fr, with_transfer = False)
            ret = Colour.__matrix_convert(fr, to, *(fr.get_params()))
            if isinstance(to, sRGB) and to.with_transfer:
                ret = (sRGB.encode_transfer(ret[0]), sRGB.encode_transfer(ret[1]), sRGB.encode_transfer(ret[2]))
            return ret

        elif isinstance(to, CIExyY):
            if isinstance(fr, sRGB) and fr.R == fr.G == fr.B == 0.:
                return (0.31272660439158, 0.32902315240275, 0.)
            else:
                if not isinstance(fr, CIEXYZ):
                    fr = CIEXYZ(fr)
                s = fr.X + fr.Y + fr.Z
                if s == 0.:
                    return (0., 0., fr.Y)
                else:
                    return (fr.X / s, fr.Y / s, fr.Y)

        elif isinstance(to, CIELAB):
            if not isinstance(fr, CIEXYZ):
                fr = CIEXYZ(fr)
            X = CIELAB.encode_transfer(fr.X / 0.95047)
            Y = CIELAB.encode_transfer(fr.Y)
            Z = CIELAB.encode_transfer(fr.Z / 1.08883)
            return (116. * Y - 16., 500. * (X - Y), 200. * (Y - Z))

        elif isinstance(to, CIE1960UCS):
            if isinstance(fr, CIEUVW):
                Y = (fr.W + 17.) / 25.
                Y *= Y * Y
                W = fr.W * 13.;
                return (fr.U / W + fr.u0, fr.V / W + fr.v0, fr.Y)
            if not isinstance(fr, CIEXYZ):
                fr = CIEXYZ(fr)
            w = fr.X + 15. * fr.Y + 3. * fr.Z
            return (4. * fr.X / w, 6. * fr.Y / w, fr.Y)

        raise Exception('Don\'t know how to convert from %s to %s' % (type(fr), type(to)))

class RGB(Colour):
    class SimpleTransferFunction(object):
        def __init__(self, gamma):
            self.gamma = gamma
            self.invgamma = 1. / gamma
        def encode(self, *rgb):
            return tuple(_pow(x, self.invgamma) for x in rgb)
        def decode(self, *rgb):
            return tuple(_pow(x, self.gamma) for x in rgb)
        def __repr__(self):
            return 'RGB.SimpleTransferFunction(%r)' % self.gamma.hex()
        def same(self, other):
            if self is other:
                return True
            if other is None or type(self) is not type(other):
                return False
            return self.gamma == other.gamma

    class RegularTransferFunction(object):
        def __init__(self, gamma, offset, slope, transition, transitioninv):
            self.gamma = gamma
            self.invgamma = 1. / gamma
            self.offset = offset
            self.slope = slope
            self.transition = transition
            self.transitioninv = transitioninv
        def encode(self, *rgb):
            def f(t):
                if t <= self.transition:
                    return self.slope * t
                return (1. + self.offset) * _pow(t, self.invgamma) - self.offset
            return tuple(f(x) for x in rgb)
        def decode(self, *rgb):
            def f(t):
                if t <= self.transitioninv:
                    return t / self.slope
                return _pow((t + self.offset) / (1. + self.offset), self.gamma)
            return tuple(f(x) for x in rgb)
        def __repr__(self):
            p = (self.gamma, self.offset, self.slope, self.transition, self.transitioninv)
            return 'RGB.RegularTransferFunction(%r, %r, %r, %r, %r)' % tuple(x.hex() for x in p)
        def same(self, other):
            if self is other:
                return True
            if other is None or type(self) is not type(other):
                return False
            if self.gamma != other.gamma:
                return False
            if self.offset != other.offset:
                return False
            if self.slope != other.slope:
                return False
            if self.transition != other.transition:
                return False
            if self.transitioninv != other.transitioninv:
                return False
            return True

    class CustomTransferFunction(object):
        def __init__(self, encoded_red, decoded_red, encoded_green = None,
                     decoded_green = None, encoded_blue = None, decoded_blue = None):
            self.encoded_red   = encoded_red
            self.decoded_red   = decoded_red
            self.encoded_green = encoded_green if encoded_green is not None else self.encoded_red
            self.decoded_green = decoded_green if decoded_green is not None else self.decoded_red
            self.encoded_blue  = encoded_blue  if encoded_blue  is not None else self.encoded_green
            self.decoded_blue  = decoded_blue  if decoded_blue  is not None else self.decoded_green
        def encode(self, R, G, B):
            return (self.encode_red(R), self.encode_green(G), self.encode_blue(G))
        def decode(self, R, G, B):
            return (self.decode_red(R), self.decode_green(G), self.decode_blue(G))
        def __repr__(self):
            p = (self.encoded_red, self.decoded_red, self.encoded_green,
                 self.decoded_green, self.encoded_blue, self.decoded_blue)
            return 'RGB.CustomTransferFunction(%r, %r, %r, %r, %r, %r)' % p

    def __init__(self, *args, with_transfer = True, transfer_function = None,
                 red = None, green = None, blue = None, white = None,
                 white_r = 1, white_g = 1, white_b = 1,
                 M = None, Minv = None, colour_space = None):
        self.with_transfer = with_transfer and transfer_function is not None
        self.transfer_function = transfer_function
        self.red   = CIExyY(red)   if red   is not None else None
        self.green = CIExyY(green) if green is not None else None
        self.blue  = CIExyY(blue)  if blue  is not None else None
        self.white = CIExyY(white) if white is not None else None
        self.white_r = white_r
        self.white_g = white_g
        self.white_b = white_b
        self.M, self.Minv = M, Minv
        self.colour_space = colour_space
        (self.R, self.G, self.B) = self.__init(args, 'RGB', True)
    def matrix_id(self):
        return list(list(r) for r in self.M)
    def get_params(self, linear = None):
        if linear is None or linear != self.with_transfer or self.transfer_function is None:
            return (self.R, self.G, self.B)
        elif linear:
            return self.transfer_function.decode(self.R, self.G, self.B)
        else:
            return self.transfer_function.encode(self.R, self.G, self.B)
    def set_params(self, R, G, B, linear = None):
        if linear is None or linear != self.with_transfer or self.transfer_function is None:
            self.R, self.G, self.B = R, G, B
        elif linear:
            self.R, self.G, self.B = return self.transfer_function.encode(R, G, B)
        else:
            self.R, self.G, self.B = return self.transfer_function.decode(R, G, B)
    def get_configuration(self):
        return {'with_transfer' : self.with_transfer,
                'transfer_function': self.transfer_function,
                'red': self.red,
                'green': self.green,
                'blue': self.blue,
                'white': self.white,
                'white_r': self.white_r,
                'white_g': self.white_g,
                'white_b': self.white_b,
                'M': self.M,
                'M': self.Minv,
                'colour_space': self.colour_space}

class sRGB(Colour):
    def __init__(self, *args, with_transfer = True):
        self.with_transfer = with_transfer
        (self.R, self.G, self.B) = self.__init(args, 'sRGB', True)
    def get_params(self, linear = None):
        if linear is None or linear != self.with_transfer:
            return (self.R, self.G, self.B)
        elif linear:
            return (sRGB.decode_transfer(self.R),
                    sRGB.decode_transfer(self.G),
                    sRGB.decode_transfer(self.B))
        else:
            return (sRGB.encode_transfer(self.R),
                    sRGB.encode_transfer(self.G),
                    sRGB.encode_transfer(self.B))
    def set_params(self, R, G, B, linear = None):
        if linear is None or linear != self.with_transfer:
            self.R, self.G, self.B = R, G, B
        elif linear:
            self.R = sRGB.encode_transfer(self.R)
            self.G = sRGB.encode_transfer(self.G)
            self.B = sRGB.encode_transfer(self.B)
        else:
            self.R = sRGB.decode_transfer(self.R)
            self.G = sRGB.decode_transfer(self.G)
            self.B = sRGB.decode_transfer(self.B)
    @staticmethod
    def encode_transfer(t):
	sign = 1
        if t < 0:
            t = -t
            sign = -1
        t = 12.92 * t if t <= 0.0031306684425217108 else.055 * t ** (1 / 2.4) - 0.055
        return t * sign
    @staticmethod
    def decode_transfer(t):
        sign = 1
        if t < 0:
            t = -t
            sign = -1
        t = t / 12.92 if t <= 0.0031306684425217108 * 12.92 else ((t + 0.055) / 1.055) ** 2.4
        return t * sign
    def get_configuration(self):
        return {'with_transfer' : self.with_transfer}

class CIExyY(Colour):
    def __init__(self, *args):
        (self.x, self.y, self.Y) = self.__init(args, 'CIExyY', False)
    def get_params(self):
        return (self.x, self.y, self.Y)
    def set_params(self, x, y, Y):
        self.x, self.y, self.Y = x, y, Y

class CIEXYZ(Colour):
    def __init__(self, *args):
        (self.X, self.Y, self.Z) = self.__init(args, 'CIEXYZ', True)
    def get_linear_params(self):
        return (self.X, self.Y, self.Z)
    def set_linear_params(self, X, Y, Z):
        self.X, self.Y, self.Z = X, Y, Z

class CIELAB(Colour):
    def __init__(self, *args):
        (self.L, self.a, self.b) = self.__init(args, 'CIELAB', False)
    def get_params(self):
        return (self.L, self.a, self.b)
    def set_params(self, L, a, b):
        self.L, self.a, self.b = L, a, b
    @staticmethod
    def encode_transfer(t):
        return t * t * t if t > 6. / 29. else (t - 4. / 29.) * 108. / 841.
    @staticmethod
    def decode_transfer(t):
        return _cbrt(t) if t > 216. / 24389. else t * 841. / 108. + 4. / 29.

class YIQ(Colour):
    def __init__(self, *args):
        (self.Y, self.I, self.Q) = self.__init(args, 'YIQ', True)
    def get_linear_params(self):
        return (self.Y, self.I, self.Q)
    def set_linear_params(self, Y, I, Q):
        self.Y, self.I, self.Q = Y, I, Q

class YDbDr(Colour):
    def __init__(self, *args):
        (self.Y, self.Db, self.Dr) = self.__init(args, 'YDbDr', True)
    def get_linear_params(self):
        return (self.Y, self.Db, self.Dr)
    def set_linear_params(self, Y, Db, Dr):
        self.Y, self.Db, self.Dr = Y, Db, Dr

class YUV(Colour):
    def __init__(self, *args):
        (self.Y, self.U, self.V) = self.__init(args, 'YUV', True)
    def get_linear_params(self):
        return (self.Y, self.U, self.V)
    def set_linear_params(self, Y, U, V):
        self.Y, self.U, self.V = Y, U, V

class YPbPr(Colour):
    def __init__(self, *args):
        (self.Y, self.Pb, self.Pr) = self.__init(args, 'YPbPr', True)
    def get_linear_params(self):
        return (self.Y, self.Pb, self.Pr)
    def set_linear_params(self, Y, Pb, Pr):
        self.Y, self.Pb, self.Pr = Y, Pb, Pr

class YCgCo(Colour):
    def __init__(self, *args):
        (self.Y, self.Cg, self.Co) = self.__init(args, 'YCgCo', True)
    def get_linear_params(self):
        return (self.Y, self.Cg, self.Co)
    def set_linear_params(self, Y, Cg, Co):
        self.Y, self.Cg, self.Co = Y, Cg, Co

class CIE1960UCS(Colour):
    def __init__(self, *args):
        (self.u, self.v, self.Y) = self.__init(args, 'CIE1960UCS', False)
    def get_params(self):
        return (self.u, self.v, self.Y)
    def set_params(self, u, v, Y):
        self.u, self.v, self.Y = u, v, Y

class CIEUVW(Colour):
    def __init__(self, *args, u0 = None, v0 = None): # TODO default to D65
        self.u0, self.v0 = u0, v0
        (self.U, self.V, self.W) = self.__init(args, 'CIEUVW', False)
    def get_params(self):
        return (self.U, self.V, self.W)
    def set_params(self, U, V, W):
        self.U, self.V, self.W = U, V, W
    def get_configuration(self):
        return {'u0' : self.u0, 'v0' : self.v0}

class CIELUV(Colour):
    def __init__(self, *args, white = None):
        self.white = CIEXYZ(white if white is not None else Illuminants[2]['D65'])
        (self.L, self.u, self.v) = self.__init(args, 'CIELUV', False)
    def get_params(self):
        return (self.L, self.u, self.v)
    def set_params(self, L, u, v):
        self.L, self.u, self.v = L, u, v
    def get_configuration(self):
        return {'white' : self.white}

class CIELChuv(Colour):
    def __init__(self, *args, white = None, one_revolution = 360.):
        self.white = CIEXYZ(white if white is not None else Illuminants[2]['D65'])
        try:
            self.one_revolution = float(one_revolution)
        except ValueError:
            self.one_revolution = float.fromhex(one_revolution)
        (self.L, self.C, self.h) = self.__init(args, 'CIELChuv', False)
    def get_params(self):
        return (self.L, self.C, self.h)
    def set_params(self, L, C, h):
        self.L, self.C, self.h = L, C, h
    def get_configuration(self):
        return {'white' : self.white, 'one_revolution' : self.one_revolution}

class YES(Colour):
    def __init__(self, *args):
        (self.Y, self.E, self.S) = self.__init(args, 'YES', True)
    def get_linear_params(self):
        return (self.Y, self.E, self.S)
    def set_linear_params(self, Y, E, S):
        self.Y, self.E, self.S = Y, E, S

Illuminants = {
    2 : {
        'A'   = CIExyY(0.447573514098910552050369915378, 0.407439444306660847328060981454, 1.),
        'B'   = CIExyY(0.348407693041403399014654951316, 0.351617234807268863594487129376, 1.),
        'C'   = CIExyY(0.310058473730255412803558101587, 0.316149707523236456196968902077, 1.),
        'D50' = CIExyY(0.345668037029273123028616510055, 0.358496838937619077825047497754, 1.),
        'D55' = CIExyY(0.332424102468830251488896010414, 0.347428039087666229445261478759, 1.),
        'D65' = CIExyY(0.312726871026564878786047074755, 0.329023206641284038376227272238, 1.),
        'D75' = CIExyY(0.299022300412497055166483050925, 0.314852737888341893679466920730, 1.),
        'E'   = CIExyY(1. / 3., 1. / 3., 1.),
        'F1'  = CIExyY(0.313062433035651010992950205036, 0.337106477918307445573731229160, 1.),
        'F2'  = CIExyY(0.372068154452825539113547392844, 0.375122558203110079144693145281, 1.),
        'F3'  = CIExyY(0.409090035308107391465171076561, 0.394117134255365986206243178458, 1.),
        'F4'  = CIExyY(0.440181095827666568620628595454, 0.403090691158138336724903183494, 1.),
        'F5'  = CIExyY(0.313756583095696484075887155996, 0.345160794752101929283583103825, 1.),
        'F6'  = CIExyY(0.377882361062687466279896852939, 0.388192885537868959122675960316, 1.),
        'F7'  = CIExyY(0.312852472915475354753311876266, 0.329174178033567632617462095368, 1.),
        'F8'  = CIExyY(0.345805753550315952971061506105, 0.358617583214377477762724311106, 1.),
        'F9'  = CIExyY(0.374105245592801061160770359493, 0.372672400924498159469067104510, 1.),
        'F10' = CIExyY(0.346086913993929323751785886998, 0.358751605952200347537939251197, 1.),
        'F11' = CIExyY(0.380537485483030235577928124258, 0.376915309293930078649026427229, 1.),
        'F12' = CIExyY(0.437023901312296902954557253906, 0.404214327891585678553809657387, 1.)
    },
    10 : {
        'A'   = CIExyY(0.451173939693730152722395132514, 0.405936604212625562482230634487, 1.),
        'B'   = CIExyY(0.349819801494100579564161535018, 0.352687989927865819250740742063, 1.),
        'C'   = CIExyY(0.310388663270034004248998371622, 0.319050711366790695766582075521, 1.),
        'D50' = CIExyY(0.347729429961154856698613002663, 0.359522508516545380441442603114, 1.),
        'D55' = CIExyY(0.334116336430253457745465084372, 0.348766090975953568786849245953, 1.),
        'D65' = CIExyY(0.313823646938709621689866935412, 0.330998985489933561510156323493, 1.),
        'D75' = CIExyY(0.299679971345752860223399238748, 0.317403239854836705102769656150, 1.),
        'E'   = CIExyY(1. / 3, 1. / 3, 1.),
        'F1'  = CIExyY(0.318098801070991199502202562144, 0.335489451474129951602520804954, 1.),
        'F2'  = CIExyY(0.379274832262508854174853922814, 0.367227934400669309145115448700, 1.),
        'F3'  = CIExyY(0.417644682102624287267644831445, 0.383124504918675723441623404142, 1.),
        'F4'  = CIExyY(0.449247699162001246087072559021, 0.390605475879083674506375700730, 1.),
        'F5'  = CIExyY(0.319739939104951298443069163113, 0.342367055369128092667807550242, 1.),
        'F6'  = CIExyY(0.386626908526034762658696308790, 0.378372201588893453116924092683, 1.),
        'F7'  = CIExyY(0.315645637312390425766039925293, 0.329508145132134222521358424274, 1.),
        'F8'  = CIExyY(0.348965563721531868424108324689, 0.359317299140994528272585739614, 1.),
        'F9'  = CIExyY(0.378258900384649654480284652891, 0.370371375730762564248976786985, 1.),
        'F10' = CIExyY(0.350893389986753234666139178444, 0.354302210111646531665030579461, 1.),
        'F11' = CIExyY(0.385435391037903751776383387551, 0.371094786781121399599214782938, 1.),
        'F12' = CIExyY(0.442654456042513022584472537346, 0.397060737666593277506166259627, 1.)
    }
}
