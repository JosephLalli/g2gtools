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

bed_fields = ["chrom", "start", "end", "strand", "motif", "annotated", "unique_reads", "multimap", "max_overhang"]
TabRecord = collections.namedtuple("TabRecord", bed_fields)

def add_TabRecords(a, b):
    assert a[0:3] == b[0:3]
    assert a[4] == b[4] # if motifs are different, then splice sites are altered and should not be merged.
    if a['annotated'] != b['annotated']:
        raise AttributeError
    annotated = max(a[5], b[5]) #if one is annotated and the other isn't, then consider the site annotated. (note: this probably shouldn't happen.)
    unique_reads = a[6] + b[6]
    multimap = (a[7] + b[7]) / 2 # ASSUMING *heuristic* that multimappers in same junction on different contigs likely mapped to both contigs, and should only be counted once
    max_overhang = max(a[8], b[8])
    return TabRecord(a[0:4] + (annotated, unique_reads, multimap, max_overhang))

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
        self.current_line_is_bed = False
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
            self.current_line_is_bed = False
            self.current_record = None
            return None

        self.current_line_is_bed = True
        elem = self.current_line.strip().split("\t")

        if not self.nitems:
            self.nitems = len(elem)
        else:
            if self.nitems != len(elem):
                raise exceptions.G2GBedError("Improperly formatted BED file")

        try:
            bed_data = {field: datum for field, datum in zip(bed_fields, elem)}

            self.current_record = TabRecord(**bed_data)
            return self.current_record
        except IndexError as ie:
            LOG.debug(ie.message)
            raise exceptions.G2GBedError("Improperly formatted BED file, line number: {0}, line: {1}".format(self.current_line_no, self.current_line))
        except ValueError as ve:
            LOG.debug(ve.message)
            raise exceptions.G2GBedError("Improperly formatted BED file, line number: {0}, line: {1}".format(self.current_line_no, self.current_line))


# TODO: kb test
def convert_tab_file(vci_file, input_file, output_file=None, reverse=False):
    """
    Convert TAB coordinates.

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

    bed_out = None
    bed_unmapped_file = None

    if output_file:
        output_file = g2g_utils.check_file(output_file, 'w')
        unmapped_file = "{0}.unmapped".format(output_file)
        bed_out = open(output_file, "w")
        bed_unmapped_file = open(unmapped_file, "w")
        LOG.info("OUTPUT FILE: {0}".format(output_file))
        LOG.info("UNMAPPED FILE: {0}".format(unmapped_file))
    else:
        input_dir, input_name = g2g_utils.get_dir_and_file(input_file)
        unmapped_file = "{0}.unmapped".format(input_name)
        unmapped_file = g2g_utils.check_file(unmapped_file, 'w')
        bed_out = sys.stdout
        bed_unmapped_file = open(unmapped_file, "w")
        LOG.info("OUTPUT FILE: stdout")
        LOG.info("UNMAPPED FILE: {0}".format(unmapped_file))

    left_right = [''] if vci_file.is_haploid() else ['_L', '_R']

    oldname_to_newname = dict()
    if reverse:
        uniq_junctions = dict()

    LOG.info("Converting BED file...")

    bed_file = BED(input_file)

    total = 0
    success = 0
    fail = 0

    for record in bed_file:

        LOG.debug("\nORIGINAL: {0}".format(str(bed_file.current_line).strip()))

        total += 1

        if total % 100000 == 0:
            LOG.info("Processed {0:,} lines".format(total))


        # If converting reference to personal, then vci and sam have the same contig names.
        # If converting personal to reference, then vci and sam have different contig names.
        # Each vci contig should be the prefix of one or more bed contigs.
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
                bed_unmapped_file.write(bed_file.current_line)
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

            elems = bed_file.current_line.rstrip().split('\t')

            LOG.debug(elems)

            elems[0] = seqid
            elems[1] = start
            elems[2] = end

            converted_junc = {field: datum for field, datum in zip(bed_fields, elems)}
            converted_junc[bed_fields[0]] = seqid
            converted_junc[bed_fields[1]] = start
            converted_junc[bed_fields[2]] = end
            converted_junc = TabRecord(**converted_junc)

            # if we're converting back to haploid, we don't want to do this twice
            junc_id = seqid+str(start)+str(end)
            if junc_id in uniq_junctions:
                uniq_junctions[junc_id] = converted_junc
            else:
                uniq_junctions[junc_id] = add_TabRecords(uniq_junctions[junc_id], converted_junc)
            if reverse:
                break
    
    for id, junc in uniq_junctions.items():
        LOG.debug("     NEW: {0}".format("\t".join(map(str, elems))))
        bed_out.write("\t".join(map(str, junc)))
        bed_out.write("\n")


    bed_out.close()
    bed_unmapped_file.close()

    LOG.info("Converted {0:,} of {1:,} records".format(success, total))
    LOG.info('TAB file converted')

