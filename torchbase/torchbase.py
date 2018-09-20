

"""
Usage:
    torchbase version [<database>] [<checkpoint>]
    torchbase run     [options] <database> <file1> [<file2>...] [--checkpoint=<checkpoint>] [--map=<mapper>]
    torchbase pull    [options] <database> [--checkpoint=<checkpoint>]
    torchbase convert [options] <fasta_file> <profile_name> <new_database_name>
    torchbase update  [options] <database_dir> <database>

Options:
    -h --help        Show this screen
    -v --verbose     Verbose logging

"""


import avro.schema
import docpie
import reference
import reference.schema
import sys
from importlib.resources import open_text




from functools import wraps



version = "0.1a"

"""command-decorator, allows subcommands to be obvious in-line, and
   binds argparse argument namespace to declared function parameters.
   ---JSP 2018-09-20
"""

command_dict = {}

schema = avro.schema.parse(open_text(reference.schema, 'package.avsc'))

def command(name):
	def command_decorator(func):
		command_dict[name] = func
		return func
	return command_decorator


@command("version")
def get_version(database=None, checkpoint=None, *a, **k):
	print(f"Torchbase v. {version}")
	if database:
		checks = reference.get_database_state(database, checkpoint)
		print(f"package {database} at checkpoint {checks[:10]}")

@command("run")
def run(file1, file2=[], *a, **k):
	pass

@command("pull")
def pull(database, *a, **k):
	pass

@command("convert")
def convert(fasta_file, profile_file, new_database_name, *a, **k):
	from convert.main import convert as c_convert
	return c_convert(fasta_file, profile_file, new_database_name)

@command("update")
def update(database_dir, database, *a, **k):
	pass

def main():
	pie = docpie.Docpie(doc=__doc__, version=None, name="torchbase", appearedonly=True)
	print(pie.preview())
	args = pie.docpie()
	args = {k.replace('>','').replace('<',''):v for k,v in args.items()}
	print(args)
	for key, command in command_dict.items():
		if args.get(key):
			exit(command_dict[key](**args))
	print(f"Torchbase v. {version}")
	print("High-performance typing schemes")
	print(__doc__)

if __name__ == '__main__':
	main()
	