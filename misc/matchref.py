#!/usr/bin/env python

import re, os, sys, codecs
from latexify import latexify

fninp = "polygly_solv.tex"
fnlib = "latex_ref_lib.bib"

class Ref:
    def __init__(self, index, author, year, title, volume = "", pages = ""):
        self.index = int(index)
        self.author = author.strip()
        s = self.author.replace(" and ", ", ")
        # replace LaTeX accent letters, such as \'e, \v{a}, and \"{\i} by a space, TODO
        s = re.sub(r"\\['`\"~\^uv](\\?[a-zA-Z]|\{\\?[a-zA-Z]\})", ", ", s)
        #s = re.sub("[^a-zA-Z, ]", ", ", s) # replace non-letter words to space
        arr = s.replace("and ", "").replace(".", ", ").split(", ")
        self.author_tokens = []
        for tok in arr:
            s = tok.strip()
            if len(s) > 1 and s != "and":
                self.author_tokens += [s,]
        self.year = int(year)
        self.title = title.strip().rstrip(".")
        # remove math from the title
        s = re.sub("\\\(.*?\\\)", "", self.title)
        s = re.sub("\$.*?\$", "", s)
        # replace non-letter symbols by a space
        # then convert to a set
        self.title_tokens = set( re.sub("[^a-zA-Z]", " ", s).split() )
        self.volume = volume.strip()
        if self.volume.isdigit():
            self.volume = int( self.volume )
        self.pages = pages.strip()
        if self.pages.isdigit():
            self.pages = int( self.pages )

    def __repr__(self):
        return "%s, %s (%s) %s %s" % (
                self.index, self.author, self.year, self.volume, self.pages)


def getref(fn):
    ''' return a list of references '''

    refstate = 0
    ls = []
    #txt = open(fn).readlines()
    txt = codecs.open(fn, encoding="utf-8").readlines()
    refstart = refend = len(txt)
    for i in range(len(txt)):
        s = txt[i]
        if refstate == 0:
            if s.startswith("\section{References}"):
                refstate = 1
                refstart = i
            continue
        else:
            ln = latexify( s.strip() )
            if ln == "": continue
            m = re.match("([0-9]+)\. .*\.", ln) # a loose pattern
            if not m: # reference ended
                refend = i
                break

            # strict pattern for volume and pages
            m = re.match("([0-9]+)\. (.*)\. ([0-9]{4})\. (.*). ([0-9a-zA-Z\- ]+): ([0-9a-zA-Z,\- ]+)\.", ln)
            if m:
                ref = Ref(m.group(1), m.group(2), m.group(3), m.group(4), m.group(5), m.group(6))
            else:
                # looser pattern for volume and pages
                m = re.match("([0-9]+)\. (.*)\. ([0-9]{4})\. (.*)\.", ln)
                if not m:
                    print "failed to get the year on %s" % ln
                    continue
                ref = Ref(m.group(1), m.group(2), m.group(3), m.group(4))
            ls += [ref,]
            #print "reference %s" % ref
    return ls, refstart, refend



def getlib(fnlib):
    ''' load the .bib file, assuming it has been latexified  '''
    txt = open(fnlib).readlines()
    n = len(txt)
    i = 0
    id = 0
    ls = []
    while i < n:
        ln = txt[i].strip()
        if ln.startswith("@article{") or ln.startswith("@book{") or ln.startswith("@incollection{"):
            tag = ln[ln.find("{")+1:-1]
            title = ""
            author = ""
            year = 0
            volume = ""
            pages = ""
            journal = ""
            while i < n:
                ln = txt[i].strip()
                if ln == "}": break
                elif ln.startswith("title = "):
                    m = re.search("title = \{(.*)\}", ln)
                    if m: title = m.group(1).replace("{", "").replace("}", "")
                elif ln.startswith("author = "):
                    m = re.search("author = \{(.*)\}", ln)
                    if m: author = m.group(1)
                elif ln.startswith("year = "):
                    m = re.search("year = \{(.*)\}", ln)
                    if m: year = m.group(1).strip()
                elif ln.startswith("volume = "):
                    m = re.search("volume = \{(.*)\}", ln)
                    if m: volume = m.group(1)
                elif ln.startswith("pages = "):
                    m = re.search("pages = \{(.*)\}", ln)
                    if m: pages = m.group(1)
                elif ln.startswith("journal = "):
                    m = re.search("journal = \{(.*)\}", ln)
                    if m: journal = m.group(1)
                i += 1
            #print id, tag, title, "by", author, "in", year, "journal", journal, "vol.", volume, "pages", pages
            #raw_input()
            ent = Ref(id, author, year, title, volume, pages)
            ent.tag = tag
            ent.journal = journal
            ls += [ent,]
            id += 1
        i += 1
    print "got %s references from %s" % (len(ls), fnlib)
    return ls



def getauthors(s):
    return s.split(", ")



def matchreflib(ref, lib):
    ''' match the references in the paper with the bibtex entries '''
    for ent in lib:
        if ent.year == ref.year and ref.volume == ent.volume and ref.pages == ent.pages:
            # check the reference authors against the library
            found = 1
            for a in ref.author_tokens:
                if a not in ent.author:
                    found = 0
                    break
            if not found: continue

            #print ent.year, "vs", ref.year, "author", ent.author
            # try to match titles
            if not ent.title_tokens.issubset(ref.title_tokens):
                print "Title mismatch for %s\n:%s\nvs\n:%s" % (ref, ent.title, ref.title)
                print ent.title_tokens, ent.volume, ent.pages
                print ref.title_tokens, ref.volume, ref.pages
                raw_input()
                continue
            return ent
    return None



def findref(ref, id):
    ''' find reference id '''
    i = id - 1
    if i >= 0 and ref[i].index == id:
        return ref[i]
    # deep search
    for i in range(len(ref)):
        if ref[i].index == id:
            return ref[i]
    return 0



def replacecite(fninp, ref, refstart, refend, fnlib):
    ''' replace citations by bibtex entries '''
    s = open(fninp).read()
    while 1:
        m = re.search("\(([0-9, \-]+)\)", s)
        if not m: break
        arr = m.group(1).split(", ")
        newarr = []
        for i in range(len(arr)):
            item = arr[i]
            p = item.find("--")
            if p >= 0: # replace range
                start = int( item[:p])
                end = int( item[p+2:] )
                newarr += range(start, end+1)
            else:
                newarr += [int(item),]

        ref_tag = []
        for item in newarr:
            r = findref(ref, item)
            if r.lib:
                ref_tag += [r.lib.tag]
        ref_tag = ", ".join(ref_tag)
        #print arr, "-->", newarr, "-->", ref_tag
        #raw_input()
        s = s[:m.start(0)] + ("\cite{%s}" % ref_tag) + s[m.end(0):]
    fnout = "ref_" + fninp
    arr = s.split("\n")
    # replace the reference section
    lib_name = os.path.splitext(fnlib)[0]
    arr = arr[:refstart] + [
            "\\bibliographystyle{abbrv}\n",
            "\\bibliography{%s}\n" % lib_name,] + arr[refend:]
    s = "\n".join(arr) + "\n"
    open(fnout, "w").write(s)
    print "saved result to %s" % fnout


if len(sys.argv) > 1:
    fninp = sys.argv[1]
if len(sys.argv) > 2:
    fnlib = sys.argv[2]

ref, refstart, refend = getref(fninp)
lib = getlib(fnlib)
for i in range(len(ref)):
    ent = matchreflib(ref[i], lib)
    if ent:
        ref[i].lib = ent
    else:
        ref[i].lib = None
        print "Ref. %s is not matched: %s" % (ref[i].index, ref[i])
replacecite(fninp, ref, refstart, refend, fnlib)


