#!/usr/bin/env python3

import re, os

fn_input_list = ["lj.html", "lj-cn.html"]
pack_dir = "pack"

class Resource:
    def __init__(self, pattern, file_index = 1, subs = []):
        self.pattern = pattern
        self.file_index = file_index
        self.subs = subs
        self.first_line_index = -1
        self.first_line = None
        self.arr = []
    
    def read_file(self, fn):
        return open(fn, encoding="utf-8").read()

    def sub_cdn(self, s, m):
        """ used to substitute resource to CDN version """
        res_fn = m.group(self.file_index)
        for sub in self.subs:
            pat = sub["pattern"]
            rep = sub["replacement"]
            m2 = re.search(pat, res_fn)
            if m2:
                s1 = s[:m.start(self.file_index)] + rep + s[m.end(self.file_index):]
                return s1
        return False

    def search(self, lines, i):
        s = lines[i]
        m = re.search(self.pattern, s)
        if m:
            # check if this file can be substituted by the CDN version
            s1 = self.sub_cdn(s, m)
            if s1:
                # replace this file by the CDN version
                lines[i] = s1
                return True

            # to combine the file
            res_fn = m.group(self.file_index)
            file = self.read_file(res_fn)
            self.arr += [file]
            lines[i] = ""

            if not self.first_line:
                self.first_line_index = i
                self.first_line = s[:m.start(self.file_index)] + "{}" + s[m.end(self.file_index):]
                #print(self.first_line)

            return True

        return False
    
    def remove_spaces(self, s):
        x = s.split("\n")
        for i in range(len(x)):
            x[i] = x[i].strip()
        return "\n".join(x)
        
    def finish(self, lines, pack_dir, filename, rmspaces = False):
        if self.first_line:
            #print(self.first_line, self.first_line_index, filename)
            lines[self.first_line_index] = self.first_line.format(filename)
        fn_output = os.path.join(pack_dir, filename)
        content = "\n".join(self.arr)
        if rmspaces: # remove redundant spaces in the combined resources
            content = self.remove_spaces(content)
        open(fn_output, "w", encoding="utf-8").write(content)
        print(f"Writing combined resource {fn_output}")


def pack(fn, pack_dir):
    basename = os.path.basename(fn)
    proj_name = basename.split(".")[0]

    proj_dir = os.path.join(pack_dir, proj_name)
    if not os.path.exists(proj_dir):
        os.makedirs(proj_dir)

    lines = open(fn, encoding="utf-8").readlines()

    js = Resource('<script .*src="(\./js.*)"></script>', 1, [
        {
            "pattern": "dygraph-combined.js",
            "replacement": "http://cdn.bootcdn.net/ajax/libs/dygraph/1.1.1/dygraph-combined.js",
        },
    ])
    css = Resource('<link .*href="(\./css/.*)"', 1)

    for i in range(len(lines)):
        if js.search(lines, i):
            continue
        if css.search(lines, i):
            continue

    js.finish(lines, proj_dir, proj_name + '_pack.js', rmspaces=False)
    css.finish(lines, proj_dir, proj_name + '_pack.css', rmspaces=False)

    fn_output = os.path.join(proj_dir, proj_name + "_pack.html")
    open(fn_output, "w", encoding="utf-8").writelines(lines)
    print(f"Writing output file {fn_output}\n")



if __name__ == "__main__":
    for fn_input in fn_input_list:
        pack(fn_input, pack_dir)
