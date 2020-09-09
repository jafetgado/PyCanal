import setuptools

with open('README.rst', 'r') as f:
    readme = f.read()

setuptools.setup(
	name='pycanal',
	version='0.1',
	author='Japheth Gado',
	author_email='japhethgado@gmail.com',
	description='Python for conservation analysis',
	long_description_content_type='text/x-rst',
	long_description=readme,
	url='https://github.com/jafetgado/pycanal',
	keywords='bioinformatics conservation-analysis,
	packages=setuptools.find_packages(),
	include_package_data=True,
	license='MIT',
	classifiers=[
		'Programming Language :: Python :: 3',
		'Operating System :: OS Independent'
				],
	install_requires=['numpy', 'pandas', 'matplotlib'],
	python_requires='>=3'
		)
