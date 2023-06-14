"""Module defining the Option class."""
import os
from pathlib import Path
from types import new_class
from typing import Callable

from lxml import etree
import xmltodict

from pyxtp.xml_editor import (
    remove_empty_text_elements,
    insert_linked_xmlfiles,
    update_xml_text,
)


class Options:
    # NOTE: if you add a class/instance attribute, please update __slots__
    __slots__ = ("_xml_tree", "_opts")

    def __init__(self, xml_file: str | Path, xml_attribs=True, set_default=False):
        """Create options object from the input xml files"""

        # parse xmlfile and replace links to the embedded xml files
        xml_file = Path(xml_file)
        self._xml_tree = insert_linked_xmlfiles(
            etree.parse(xml_file), base_path=f"{xml_file.parent}"
        )
        dict_data = xmltodict.parse(
            etree.tostring(self._xml_tree), xml_attribs=xml_attribs
        )
        name = xml_file.stem
        self._opts = getattr(make_options(name, dict_data, set_default).options, name)

    def __getattr__(self, name: str):
        if name in self.__slots__:
            return getattr(self, name)
        return getattr(self._opts, name)

    def __setattr__(self, name: str, value):
        if name in self.__slots__:
            super().__setattr__(name, value)
        else:
            setattr(self._opts, name, value)

    def __dir__(self):
        from_opts = [attr for attr in dir(self._opts) if not attr.startswith("_")]
        return [*self.__slots__, *from_opts]

    def _write_xml(self, filename: str) -> None:
        """update and clean the xml and then write the input file

        Args:
            filename (str): filename in the
        """

        # create a dict of the user providedoptions
        options = self._opts.to_flat_dict()

        # update the options from the user provided input
        xml_tree_cpy = update_xml_text(self._xml_tree, options)

        # remove the empty elements
        xml_tree_cpy = remove_empty_text_elements(xml_tree_cpy)

        # write the file
        xml_tree_cpy.write(filename)


class _Opts_t:
    __slots__ = ()
    def to_dict(self) -> dict:
        """Convert a config object to a dictionary, skip keys that are not set"""
        res = {}
        for attr in dir(self):
            if attr.startswith("_"):
                continue
            value = getattr(self, attr)
            if value is None or isinstance(value, Callable):
                continue
            if isinstance(value, _Opts_t):  # recurse
                if _val := value.to_dict():  # drop empty
                    res[attr] = _val
            else:  # leaf node
                res[attr] = value
        return res

    def to_flat_dict(self) -> dict:
        """Transform options object into a flat dictionary.

        Returns:
            dict: dictionary maps the namespace
        """

        def _flatten_dict(data: dict, path: str = "") -> dict:
            """Flatten a nested dictionary.

            Args:
                data (dict): input dictionary
                path (str, optional): basepath in the structure. Defaults to ''.

            Returns:
                dict: dictionary maps the namespace
            """
            output = dict()
            for key, value in data.items():
                if isinstance(value, dict):
                    output.update(_flatten_dict(value, path=os.path.join(path, key)))
                else:
                    if not key.startswith("_"):
                        output[os.path.join(path, key)] = value
            return output

        return _flatten_dict(self.to_dict())


def make_options(name: str, xml_dict: dict, set_default: bool = True):
    """Create a config object out of an XML dictionary

    Parameters
    ----------
    name : str
      Name for options type

    xml_dict : dict
      Dictionary generated from the XML configuration files

    set_default : bool (default: True)
      Whether to set default values in the config object

    Returns
    -------
    name_t

    """
    # detect leaf, return early w/ default
    children = [k for k in xml_dict if not k.startswith("@")]
    if len(children) == 0:  # leaf node
        default = xml_dict.get("@default", None)
        if default == "OPTIONAL":
            return None
        if not set_default:
            return ''
        return default

    # create type
    attrs, slots = {}, {}
    for key, value in xml_dict.items():
        if isinstance(value, dict):
            doc = [v for k, v in value.items() if k.startswith("@")]
            slots[key] = "\n".join(doc)
            attrs[key] = make_options(key, value, set_default)
    cls = new_class(
        f"{name}_t",
        bases=(_Opts_t,),
        exec_body=lambda ns: ns.update({"__slots__": slots}),
    )
    obj = cls()

    # assign slots & defaults
    for k, v in attrs.items():  # no @keys here
        if isinstance(v, _Opts_t) or set_default:
            setattr(obj, k, v)
        elif not set_default:
            setattr(obj, k, '')
        else:
            raise ValueError('Something went wrong')
    return obj
