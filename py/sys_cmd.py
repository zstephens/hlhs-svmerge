import os

def exe(cmd):
	print cmd
	os.system(cmd)

def exists_and_is_nonZero(fn):
	if os.path.isfile(fn):
		if os.path.getsize(fn) > 0:
			return True
	return False

def dir_format(d):
	if d[-1] != '/':
		return d+'/'
	else:
		return d
		