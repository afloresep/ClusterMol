from setuptools import setup, find_packages

setup(
    name='molclust',
    version='0.1.0',
    description='Molecular clustering and UMAP-based analysis tool',
    author='Alejandro Flores',
    packages=find_packages(),
    install_requires=[
        'pandas',
        'scikit-learn',
        'umap-learn',
        'openpyxl',
        'rdkit',
    ],
    entry_points={
        'console_scripts': [
            'molclust=molclust.main:main', # For CLI entry point 
        ],
    },
    include_package_data=True,
    zip_safe=False
)