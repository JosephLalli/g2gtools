# -*- coding: utf-8 -*-

#
# Collection of functions related to BED files
#
# 0-based
#

import collections
import sys

from . import compat
from . import exceptions
from . import g2g
from . import g2g_utils
from . import vci


junc_fields = ["chrom", "start", "end", "name", "score", "strand", "extra"]
JuncRecord = collections.namedtuple("JuncRecord", junc_fields)

def add_JuncRecords(a, b):
    assert a[0:3] == b[0:3]
    both = str(int(a[4]) + int(b[4]))
    return JuncRecord(a.chrom, a.start, a.end, a.name, both, a.strand, a.extra)


LOG = g2g.get_logger()


class BED(object):
    """
    Simple BED object for parsing and iterating BED files.

    Supports transparent gzip decompression.
    """
    def __init__(self, filename):
        if not filename:
            raise exceptions.G2GBedError("A filename must be supplied")

        self.filename = filename
        self.current_line = None
        self.current_line_is_junc = False
        self.current_record = None
        self.reader = g2g_utils.open_resource(filename)
        self.nitems = None
        self.current_line_no = 0

    def __iter__(self):
        return self

    def __next__(self):

        if compat.is_py2:
            self.current_line = self.reader.next()
            self.current_line_no += 1

            while self.current_line and len(self.current_line.strip()) == 0:
                self.current_line = self.reader.next()
                self.current_line_no += 1

        else:
            self.current_line = g2g_utils.s(self.reader.__next__())
            self.current_line_no += 1

            while self.current_line and len(self.current_line.strip()) == 0:
                self.current_line = self.reader.__next__()
                self.current_line_no += 1

        if self.current_line.startswith("track"):
            self.current_line = self.current_line.strip()
            self.current_line_is_junc = False
            self.current_record = None
            return None

        self.current_line_is_junc = True
        elem = self.current_line.strip().split("\t")

        if not self.nitems:
            self.nitems = len(elem)
        else:
            if self.nitems != len(elem):
                raise exceptions.G2GBedError("Improperly formatted BED file")

        try:
            junc_data = {'chrom': elem[0],
                         'start': int(elem[1]),
                         'end': int(elem[2]),
                         'name': elem[3] if self.nitems > 3 else None,
                         'score': elem[4] if self.nitems > 4 else None,
                         'strand':  elem[5] if self.nitems > 5 else None,
                         'extra': '\t'.join(elem[6:]) if self.nitems > 6 else None}

            self.current_record = JuncRecord(**junc_data)
            return self.current_record
        except IndexError as ie:
            LOG.debug(ie.message)
            raise exceptions.G2GBedError("Improperly formatted BED file, line number: {0}, line: {1}".format(self.current_line_no, self.current_line))
        except ValueError as ve:
            LOG.debug(ve.message)
            raise exceptions.G2GBedError("Improperly formatted BED file, line number: {0}, line: {1}".format(self.current_line_no, self.current_line))


# TODO: kb test
def convert_junc_file(vci_file, input_file, output_file=None, reverse=False):
    """
    Convert JUNC coordinates.

    :param vci_file: VCI input file
    :type vci_file: :class:`.chain.ChainFile`
    :param input_file: the input BED file
    :type input_file: string
    :param output_file: the output BED file
    :type output_file: string
    :param reverse: reverse direction of original file
    :type reverse: boolean
    :return:
    """

    if isinstance(vci_file, vci.VCIFile):
        LOG.info("VCI FILE: {0}".format(vci_file.filename))
        LOG.info("VCI FILE IS DIPLOID: {0}".format(vci_file.is_diploid()))
    else:
        vci_file = g2g_utils.check_file(vci_file)
        vci_file = vci.VCIFile(vci_file)
        LOG.info("VCI FILE: {0}".format(vci_file.filename))
        LOG.info("VCI FILE IS DIPLOID: {0}".format(vci_file.is_diploid()))
        vci_file.parse(reverse)

    input_file = g2g_utils.check_file(input_file)
    LOG.info("INPUT FILE: {0}".format(input_file))

    junc_out = None
    junc_unmapped_file = None

    if output_file:
        output_file = g2g_utils.check_file(output_file, 'w')
        unmapped_file = "{0}.unmapped".format(output_file)
        junc_out = open(output_file, "w")
        junc_unmapped_file = open(unmapped_file, "w")
        LOG.info("OUTPUT FILE: {0}".format(output_file))
        LOG.info("UNMAPPED FILE: {0}".format(unmapped_file))
    else:
        input_dir, input_name = g2g_utils.get_dir_and_file(input_file)
        unmapped_file = "{0}.unmapped".format(input_name)
        unmapped_file = g2g_utils.check_file(unmapped_file, 'w')
        junc_out = sys.stdout
        junc_unmapped_file = open(unmapped_file, "w")
        LOG.info("OUTPUT FILE: stdout")
        LOG.info("UNMAPPED FILE: {0}".format(unmapped_file))

    left_right = [''] if vci_file.is_haploid() else ['_L', '_R']

    oldname_to_newname = dict()
    if reverse:
        uniq_junctions = dict()

    LOG.info("Converting BED file...")

    junc_file = BED(input_file)

    total = 0
    success = 0
    fail = 0

    for record in junc_file:

        LOG.debug("\nORIGINAL: {0}".format(str(junc_file.current_line).strip()))

        total += 1

        if total % 100000 == 0:
            LOG.info("Processed {0:,} lines".format(total))


        # If converting reference to personal, then vci and sam have the same contig names.
        # If converting personal to reference, then vci and sam have different contig names.
        # Each vci contig should be the prefix of one or more junc contigs.
        # Name to id converts input sam contig names to output tids. 
        # So, both contig_L and contig_R should be assigned to the same id, contig's id

        if vci_file.is_diploid() and reverse:
            if record.chrom not in oldname_to_newname.keys():
                matching_contigs = [c for c in vci_file.contigs if record.chrom[:len(c)] == c]
                if len(matching_contigs) == 1:
                    oldname_to_newname[record.chrom] = matching_contigs[0]
                else:
                    raise AttributeError


        for lr in left_right:
            if reverse:
                seqid = record.chrom
            else:
                seqid = "{}{}".format(record.chrom, lr)
            
            mappings = vci_file.find_mappings(seqid, record.start - 1, record.end)

            # unmapped
            if mappings is None:
                # LOG.info((seqid, record.start - 1, record.end),"\tFail due to no mappings")
                # LOG.info (vci_file.contigs)
                junc_unmapped_file.write(junc_file.current_line)
                fail += 0
                continue
            else:
                LOG.debug("{0} mappings found".format(len(mappings)))
                if reverse:
                    seqid = oldname_to_newname[record.chrom]

            success += 1
            start = mappings[0].to_start + 1
            end = mappings[-1].to_end

            LOG.debug("({0}, {1}) => ({2}, {3})".format(record.start - 1, record.end, start, end))

            elems = junc_file.current_line.rstrip().split('\t')

            LOG.debug(elems)

            elems[0] = seqid
            elems[1] = start
            elems[2] = end

            converted_junc = {field: datum for field, datum in zip(junc_fields, elems)}
            converted_junc[junc_fields[0]] = seqid
            converted_junc[junc_fields[1]] = start
            converted_junc[junc_fields[2]] = end
            converted_junc = JuncRecord(**converted_junc)

            # if we're converting back to haploid, we don't want to do this twice
            junc_id = seqid+str(start)+str(end)
            if junc_id not in uniq_junctions:
                uniq_junctions[junc_id] = converted_junc
            else:
                uniq_junctions[junc_id] = add_JuncRecords(uniq_junctions[junc_id], converted_junc)
            if reverse:
                break
    
    for id, junc in uniq_junctions.items():
        LOG.debug("     NEW: {0}".format("\t".join(map(str, elems))))
        junc_out.write("\t".join(map(str, junc)))
        junc_out.write("\n")


    junc_out.close()
    junc_unmapped_file.close()

    LOG.info("Converted {0:,} of {1:,} records".format(success, total))
    LOG.info('JUNC file converted')