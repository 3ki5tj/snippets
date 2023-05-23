#!/usr/bin/env python3


class CHeader:
  def __init__(self, fn):
    self.fn = fn
    self.raw = open(fn).readlines()
    self.s = self.trim_lines(self.raw)
    self.include_statement = '#include "%s"' % fn
    self.included = False


  def trim_lines(self, lines):
    s = lines[:]
    
    for start in range(5):
      if (  s[start].startswith('#ifndef ')
        and s[start+1].startswith('#define ') ):
        s[start] = ''
        s[start+1] = ''
        break
    
    for end in range(len(s)-1, 0, -1):
      if s[end].startswith('#endif '):
        s[end] = ''

    return s


  def replace_include(self, another_header):

    i = 0
    while i < len(self.s):
      ln = self.s[i]
      if another_header.include_statement in ln:
        if another_header.included:
          self.s[i] = ''
        else:
          self.s = self.s[:i] + another_header.s + self.s[i+1:]
          another_header.included = True
      i += 1


  def get_macro_name(self, fn):

    p = fn.rfind('/')
    if p >= 0:
      fn = fn[p+1:]

    s = fn.replace(".h", "_H__").upper()

    return s


  def save(self, fn, with_macro = True):

    macro_name = self.get_macro_name(fn)

    if with_macro:
      s = [
        '#ifndef %s\n' % macro_name,
        '#define %s\n' % macro_name,
      ]
      s += self.s
      s += [
        '#endif /* defined(%s) */\n' % macro_name,
      ]
    else:
      s = self.s

    open(fn, "w").writelines(s)



def pack_rng():
  mt = CHeader('rng_engine_mt.h')
  pcg = CHeader('rng_engine_pcg.h')
  mgr = CHeader('rng_engine_manager.h')
  rng = CHeader('rng.h')
  mgr.replace_include(mt)
  mgr.replace_include(pcg)
  rng.replace_include(mgr)
  rng.save('rng_pack.h')

pack_rng()
