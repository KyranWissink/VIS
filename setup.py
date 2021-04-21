from setuptools import setup
def get_contents(*args):
    """Get the contents of a file relative to the source distribution directory."""
    with codecs.open(get_absolute_path(*args), 'r', 'UTF-8') as handle:
        return handle.read()

def get_requirements(*args):
    """Get requirements from pip requirement files."""
    requirements = set()
    with open(get_absolute_path(*args)) as handle:
        for line in handle:
            # Strip comments.
            line = re.sub(r'^#.*|\s#.*', '', line)
            # Ignore empty lines
            if line and not line.isspace():
                requirements.add(re.sub(r'\s+', '', line))
    return sorted(requirements)

setup(
    name='VIS',
    version='0.1.0',
    install_requires=get_requirements('requirements.txt'),
    url='https://github.com/KyranWissink/VIS',
    license='',
    author='Kyran Wissink',
    author_email='',
    description='VIS: HGVS variant interpretation using SPLICEAI',
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "Operating System :: OS Independent",
    ],
    long_description=get_contents('README.md'),
    python_requires=">=3.7",
)