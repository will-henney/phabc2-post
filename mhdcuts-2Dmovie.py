import subprocess, sys, os

try:
    runid = sys.argv[1]
    imageprefix = "hsv-xtd-bbb-cuts-%s" % (runid)
except:
    print "Usage: %s RUNID" % (os.path.basename(sys.argv[0]))
encode_exec = "mencoder"
encode_args = "-ovc lavc -lavcopts vbitrate=5000:vcodec=wmv2 -mf type=png:fps=15"
file_args = r"-o %s.avi mf://%s-\?\?\?\?.png" % (imageprefix, imageprefix)
cmd = "%s %s %s" % (encode_exec, encode_args, file_args)
subprocess.Popen(cmd, shell=True).wait()
print "Written %s.avi" % (imageprefix)

