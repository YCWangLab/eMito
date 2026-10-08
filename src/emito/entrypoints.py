"""Installed entry points for independent modules."""
import sys
from .cli import main

def prepare():
    return main(["prepare", *sys.argv[1:]])

def taxa_generate():
    return main(["taxa-generate", *sys.argv[1:]])

def node_generate():
    return main(["node-generate", *sys.argv[1:]])

def access():
    return main(["access", *sys.argv[1:]])

def collapse():
    return main(["collapse", *sys.argv[1:]])

def merge():
    return main(["merge", *sys.argv[1:]])
