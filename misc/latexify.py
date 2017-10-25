#!/usr/bin/env python

import codecs, os, sys, re

# converting .bib file to LaTeX format
# may work for other files as well


#https://www.johndcook.com/blog/greek_letters/
chmap = {
  u"\u0391": r"$\Alpha$",
  u"\u0392": r"$\Beta$",
  u"\u0393": r"$\Gamma$",
  u"\u0394": r"$\Delta$",
  u"\u0395": r"E",
  u"\u0396": r"Z",
  u"\u0397": r"H",
  u"\u0398": r"$\Theta$",
  u"\u0399": r"I",
  u"\u039A": r"K",
  u"\u039B": r"$\Lambda$",
  u"\u039C": r"M",
  u"\u039D": r"N",
  u"\u039E": r"$\Xi$",
  u"\u039F": r"O",
  u"\u03A0": r"$\Pi$",
  u"\u03A1": r"P",
  u"\u03A3": r"$\Sigma$",
  u"\u03A4": r"$\Tau$",
  u"\u03A5": r"$\Upsilon$",
  u"\u03A6": r"$\Phi$",
  u"\u03A7": r"X",
  u"\u03A8": r"$\Psi$",
  u"\u03A9": r"$\Omega$",
  #
  u"\u03b1": r"$\alpha$",
  u"\u03b2": r"$\beta$",
  u"\u03b3": r"$\gamma$",
  u"\u03b4": r"$\delta$",
  u"\u03b5": r"$\epsilon$",
  u"\u03b6": r"$\zeta$",
  u"\u03b7": r"$\eta$",
  u"\u03b8": r"$\theta$",
  u"\u03b9": r"$\iota$",
  u"\u03ba": r"$\kappa$",
  u"\u03bb": r"$\lambda$",
  u"\u03bc": r"$\mu$",
  u"\u03bd": r"$\nu$",
  u"\u03be": r"$\xi$",
  u"\u03bf": r"o",
  u"\u03c0": r"$\pi$",
  u"\u03c1": r"$\rho$",
  u"\u03c2": r"$\varsigma$",
  u"\u03c3": r"$\sigma$",
  u"\u03c4": r"$\tau$",
  u"\u03c5": r"$\upsilon$",
  u"\u03c6": r"$\varphi$",
  u"\u03c7": r"$\chi$",
  u"\u03c8": r"$\psi$",
  u"\u03c9": r"$\omega$",
  u"\u03d5": r"$\phi$",
  # https://tex.stackexchange.com/tags/accents/info
  u"\u00b0": r"$^{\circ}$",
  u"\u00b7": r"$\cdot$",
  u"\u00d7": r"$\times$",
  u"\u00e1": r"\'{a}",
  u"\u00e5": r"\aa",
  u"\u00e8": r"\`{e}",
  u"\u00e9": r"\'{e}",
  u"\u00eb": r'\"{e}',
  u"\u00ed": r"\'{i}",
  u"\u00ec": r"\`{i}",
  u"\u00ef": r'\"{i}',
  u"\u00f3": r"\'{o}",
  u"\u0107": r"\'{c}",
  u"\u00e3": r'\~{a}',
  u"\u00e4": r'\"{a}',
  u"\u00e7": r'\c{c}',
  u"\u00ee": r'\^{i}',
  u"\u00f6": r'\"{o}',
  u"\u00f8": r'{\o}',
  u"\u00fa": r"\'u",
  u"\u00fc": r'\"u',
  u"\u010d": r'\v{c}',
  u"\u011b": r'\v{e}',
  u"\u0131": r'{\i}',
  u"\u0142": r'{\l}',
  u"\u0144": r"\'{n}",
  u"\u014d": r"\={o}",
  u"\u0161": r'\v{s}',
  u"\u025b": r"$\varepsilon",
  u"\u2212": r"--",
  u"\u2009": r"\;",
  u"\u2010": r"-",
  u"\u2013": r"-",
  u"\u2014": r"--",
  u"\u2018": r"`",
  u"\u2019": r"'",
  u"\u2113": r"$\ell$",
  u"\u2022": r"$\cdot$",
  u"\u2218": r"$\circ$",
  u"\u22c5": r"$\cdot$",
  u"\u221e": r"$\infty$",
  u"\u22ef": r"\dots",
  u"\u201c": r"``",
  u"\u201d": r"''",
  u"\u00c5": r"\AA",
  u"\u00a9": r"\copyright",
  u"\u2122": r"\texttrademark",
  u"\u21cc": r"$\rightleftharpoons$",
  u"\u00b1": r"$\pm",
  u"\u2a7d": r"$\le$", # r"$\leqslant$",
  u"\u2a7e": r"$\ge$", # r"$\geqslant$",
  u"\u2264": r"$\le$",
  u"\u2265": r"$\ge$",
  u"\u221d": r"$\propto$",
  u"\u223c": r"$\sim$",
  u"\u2248": r"$\approx$",
  u"\u2020": r"$\dagger$",
  u"\u2032": r"$'$",
  u"\u3008": r"$\langle$",
  u"\u3009": r"$\rangle$",
  u"\u2192": r"$\to$",
}

combine_accents = {
  u"\u0300": r"\`",
  u"\u0301": r"\'",
  u"\u0302": r"\^",
  u"\u0303": r"\~",
  u"\u0304": r"\=",
  u"\u0305": r"\=",
  u"\u0306": r"\u",
  u"\u0307": r"\.",
  u"\u0308": r'\"'
}

misc_repl = {
  r"\&amp;": r"\&",
  r"~": "$\sim$",
}

def latexify(s):
    ''' convert fninp to LaTeX format '''

    # handle combining accents
    for c in combine_accents:
        mc = combine_accents[c]
        while 1:
            p = s.find(c)
            if p >= 0:
                print s[p+1]+s[p], "-->", mc + "{" + s[p+1] + "} at position", p
                s = s[:p] + mc + "{" + s[p+1] + "}" + s[p+2:]
            else:
                break

    # normal substitution
    for c in chmap:
        s = s.replace(c, chmap[c])

    for c in misc_repl:
        s = s.replace(c, misc_repl[c])
    return s

    ## use the following loop to detect strange strings
    # unknown = []
    #for c in s:
    #    if ord(c) >= 128 and c not in chmap:
    #        if c not in unknown:
    #            unknown += [c,]
    #            print "%d: [%s] in [%04x]" % (len(unknown), c, ord(c))


def rmlanguage(s):
    ''' remove language entries in .bib file '''
    s = re.sub(r"\n\s+?language = {(en|eng)},\s*?\n", "\n", s)
    return s
    

fninp = "ref_lib.bib"
if len(sys.argv) > 1:
    fninp = sys.argv[1]
fnout = "latex_" + fninp

s = codecs.open(fninp, encoding="utf-8").read()
s = latexify(s)
s = rmlanguage(s)
codecs.open(fnout, mode="w", encoding="utf-8").write(s)
print "save file to %s" % fnout
