from setuptools import setup, find_packages
from setuptools.extension import Extension

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name='sideEffects',
    description='sideEffects',
    py_modules=['sideEffects', 'accel'],
    install_requires=requirements,
    zip_safe=False,
    entry_points={
        'console_scripts':[
            'sideEffects = sideEffects.processing:main'
        ]
    },
    ext_modules=[
        Extension("sideEffects.accel", sources = ["sideEffects/extensions/accel.cpp"])
    ],
    test_suite="tests"
)
