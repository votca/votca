from functools import reduce
import os, os.path

import pytest
from lxml import etree
import xmltodict

from pyxtp.options import Options
from pyxtp.xml_editor import insert_linked_xmlfiles


@pytest.fixture(scope="module")
def dftgwbse_xml():
    votcashare = os.environ.get("VOTCASHARE")
    if votcashare is None:
        msg = (
            "pyxtp: cannot find Votca installation, "
            "please set the VOTCASHARE environment variable"
        )
        raise RuntimeError(msg)
    return f"{votcashare}/xtp/xml/dftgwbse.xml"


class TestOptions:
    @pytest.fixture
    def opts(self, dftgwbse_xml):
        return Options(dftgwbse_xml, set_default=True)

    @pytest.fixture
    def xml_dict(self, dftgwbse_xml):
        path = dftgwbse_xml
        parent = os.path.dirname(path)
        xml_tree = insert_linked_xmlfiles(etree.parse(path), base_path=parent)
        return xmltodict.parse(etree.tostring(xml_tree), xml_attribs=True)

    def test_completion(self, opts, xml_dict):
        # check tab completion of properties
        expected = set(attr for attr in dir(opts) if not attr.startswith("to_"))
        from_xml = list(
            k for k in xml_dict["options"]["dftgwbse"] if not k.startswith("@")
        )
        if expected != {*opts.__slots__, *from_xml}:
            raise AssertionError("Error in test_completion")

    @pytest.mark.parametrize(
        "attr, val", [("job_name", "options_test"), ("mpsfile", "somefile")]
    )
    def test_attr_valid(self, opts, attr, val):
        setattr(opts, attr, val)

    def test_attr_invalid(self, opts):
        with pytest.raises(AttributeError, match="jobname"):
            opts.jobname = "options_test"

    @pytest.mark.parametrize(
        "attr, expected", [("job_name", "system"), (("dftpackage", "spin"), "1")]
    )
    def test_attr_default(self, opts, attr, expected):
        if isinstance(attr, str):
            val = getattr(opts, attr)
        else:
            val = reduce(getattr, attr, opts)
        if expected != val:
            raise AssertionError("Error in test_attr_default")