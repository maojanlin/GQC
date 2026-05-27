import sys
import logging
import pybedtools
from GQC import bedtoolslib

logger = logging.getLogger(__name__)

# reminders: in aligndata, (1) all coordinates are 1-based, (2) strand is "+" or "-", (3) querystart is the query's lower coordinate,
# so doesn't correspond to targetstart if alignment is on the reverse strand

def add_structural_error(structural_errors:list, chrom:str, start:int, end:int, errortype:str, query:str, query1, query2, widestart, wideend, netdiff, strand:str):

    structural_errors.append({
        "chrom": chrom,
        "start": start,
        "end": end,
        "line": chrom + "\t" + str(start) + "\t" + str(end) + "\t" + errortype + "\t" + query + "\t" + str(query1) + "\t" + str(query2) + "\t" + str(widestart) + "\t" + str(wideend) + "\t" + str(netdiff) + "\t" + strand + "\n"
    })

    return 0

def write_structural_errors(aligndata:list, refobj, queryobj, outputdict, bmstats, args, excludedbedobj=None)->str:

    aligndict = {}
    current_align = None
    structural_errors = []
    for align in sorted(aligndata, key=lambda a: (a["target"], a["targetstart"], a["targetend"])):
        refentry = align["target"]
        if refentry not in aligndict:
            aligndict[refentry] = [align]
        else:
            aligndict[refentry].append(align)
        query = align["query"]
        refstart = align["targetstart"]
        refend = align["targetend"]
        querystart = align["querystart"]
        queryend = align["queryend"]
        strand = align["strand"]
        if current_align is not None:
            refdiff = refstart - current_align["targetend"]
            if refentry == current_align["target"] and query == current_align["query"] and strand == current_align["strand"]:
                if strand == "+":
                    querydiff = querystart - current_align["queryend"]
                    query1 = current_align["queryend"]
                    query2 = querystart
                else:
                    querydiff = queryend - current_align["querystart"]
                    query1 = querystart
                    query2 = current_align["queryend"]

                netdiff = querydiff - refdiff
                if refdiff < querydiff: # refdiff less than querydiff (insertion), netshift positive
                    if refdiff > 0:
                        add_structural_error(structural_errors, refentry, current_align["targetend"] - 1, refstart, "SameContigInsertion", query, query1, query2, current_align["targetend"], refstart, netdiff, strand)
                    else:
                        add_structural_error(structural_errors, refentry, refstart - 1, current_align["targetend"], "SameContigInsertion", query, query1, query2, current_align["targetend"], refstart, netdiff, strand)
                else: # refdiff greater than than querydiff (deletion), netshift negative
                    if refdiff > 0:
                        add_structural_error(structural_errors, refentry, current_align["targetend"] - 1, refstart, "SameContigDeletion", query, query1, query2, current_align["targetend"], refstart, netdiff, strand)
                    else:
                        add_structural_error(structural_errors, refentry, refstart - 1, current_align["targetend"], "SameContigDeletion", query, query1, query2, current_align["targetend"], refstart, netdiff, strand)

            elif refentry == current_align["target"]: # strand switch or new contig:
                queryentries = query + "/" + current_align["query"]
                strands = strand + "/" + current_align["strand"]
                if refdiff > 0:
                    add_structural_error(structural_errors, refentry, current_align["targetend"] - 1, refstart, "BetweenContigDeletion", queryentries, ".", ".", current_align["targetend"], refstart, "NA", strands)
                else:
                    add_structural_error(structural_errors, refentry, refstart - 1, current_align["targetend"], "BetweenContigInsertion", queryentries, ".", ".", current_align["targetend"], refstart, "NA", strands)
        current_align = align

    excluded_indices = set()
    if structural_errors and excludedbedobj is not None and len(excludedbedobj) > 0:
        structvarbedstring = ""
        for index, structural_error in enumerate(structural_errors):
            bed_start, bed_end = bedtoolslib.normalize_bed_interval(structural_error["start"], structural_error["end"])
            structvarbedstring = structvarbedstring + structural_error["chrom"] + "\t" + str(bed_start) + "\t" + str(bed_end) + "\t" + str(index) + "\n"
        structvarbedobj = pybedtools.BedTool(structvarbedstring, from_string=True)
        excludedsvs = structvarbedobj.intersect(excludedbedobj, wa=True)
        for excludedsv in excludedsvs:
            excluded_indices.add(int(excludedsv.name))
        logger.info("Excluded " + str(len(excluded_indices)) + " structural variants overlapping excluded benchmark regions")

    excludedfh = None
    if "excludedstructvariantbed" in outputdict.keys():
        excludedfh = open(outputdict["excludedstructvariantbed"], "w")
    try:
        with open(outputdict["structvariantbed"], "w") as sfh:
            for index, structural_error in enumerate(structural_errors):
                if index in excluded_indices:
                    if excludedfh:
                        excludedfh.write(structural_error["line"])
                else:
                    sfh.write(structural_error["line"])
    finally:
        if excludedfh:
            excludedfh.close()

    return 0

