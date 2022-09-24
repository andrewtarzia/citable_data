#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to write meta data for each folder.

Author: Andrew Tarzia

"""

import logging
import json
import pathlib
import os


class Choice:
    def __init__(self, string):
        if string.lower() == "y":
            self._result = True
        elif string.lower() == "n":
            self._result = False
        else:
            logging.info("you did not provide y/n - try again:")
            self._result = Choice(input("try again, y/n"))

    def get_choice(self):
        return self._result

    def __str__(self) -> str:
        return f"{self.__class__.__name__}({self._result})"

    def __repr__(self) -> str:
        return str(self)


class FileChoice(Choice):
    def __init__(self):
        string = input(
            "tell me. Is it a structure file (s), "
            "input file (i), output file (o), "
            "data file (d), figure (f), run file (r), "
            "or should it be ignored (n)?"
        )

        if string.lower() == "s":
            self._result = "structure"
        elif string.lower() == "o":
            self._result = "output"
        elif string.lower() == "d":
            self._result = "data"
        elif string.lower() == "i":
            self._result = "input"
        elif string.lower() == "f":
            self._result = "figure"
        elif string.lower() == "r":
            self._result = "runfile"
        elif string.lower() == "n":
            self._result = "ignore"
        else:
            logging.info(
                "you did not provide any of s/o/i/d/f/r/n - try again:"
            )
            self._result = Choice(input("try again, s/o/i/d/f/r/n"))


def ignored_subdirs():
    return ("/data/atarzia/projects/citable_data", ".git")


def get_mapping_rules():
    dir_path = pathlib.Path(os.path.dirname(os.path.realpath(__file__)))
    with open(dir_path / "mapping_rules.json", "r") as f:
        maps = json.load(f)
    return maps


def update_mapping_rules(new_key, new_val):
    maps = get_mapping_rules()
    maps[new_key] = new_val

    dir_path = pathlib.Path(os.path.dirname(os.path.realpath(__file__)))
    with open(dir_path / "mapping_rules.json", "w") as f:
        json.dump(maps, f)

    return get_mapping_rules()


def main():

    dir_path = pathlib.Path(os.path.dirname(os.path.realpath(__file__)))

    entries = dir_path.iterdir()
    entries = [
        i
        for i in sorted(entries, key=lambda entry: entry.is_file())
        if i.name not in ignored_subdirs() and i.is_dir()
    ]
    entries_count = len(entries)
    logging.info(f"there are {entries_count} projects")

    mappings = get_mapping_rules()
    for entry in entries:
        json_file = entry / f"{entry.name}_meta.json"
        if os.path.exists(json_file):
            continue
        logging.info(f"filling {entry} as {entry.name}")

        contents = {
            "structure": [],
            "input": [],
            "output": [],
            "data": [],
            "figure": [],
            "runfile": [],
        }

        # I want to add to contents, all input, output and data
        # files.
        for subdir, dirs, files in os.walk(entry):
            for fi in files:
                path = pathlib.Path(os.path.join(subdir, fi))

                # Handle unknown file extensions.
                try:
                    file_type = mappings[path.suffix]
                except KeyError:
                    logging.info(
                        f"I do not know what a {path.suffix} is"
                    )
                    choice = FileChoice().get_choice()
                    mappings = update_mapping_rules(
                        new_key=path.suffix,
                        new_val=choice,
                    )
                    file_type = mappings[path.suffix]

                if file_type in ("ignore"):
                    continue

                contents[file_type].append(
                    str(path).replace(str(entry), "")
                )

        with open(json_file, "w") as f:
            json.dump(contents, f, indent=4)


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s | %(levelname)s | %(message)s",
    )
    main()
