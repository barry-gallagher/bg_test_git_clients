import os
import re
import logging

LOGGER = logging.getLogger('proc.enc_sup')


sorind_list = ["US,US,graph,DD1", "US,US,graph,graph,DD1", "US,US,reprt,DD2", "US,US,US,graph,DD1"]
ident_breakdown = {'g': 'graph', 'r': 'reprt'}
format_string = f"US,U[1-4S],(?:{('|').join(ident_breakdown.values())}),"
ident_chars = f"[{('').join(ident_breakdown.keys())}][a-z]" + '{4}'
naming_standard = re.compile(format_string)  # , re.IGNORECASE)
index_breakdown = {
    0: (2, re.compile("US", re.IGNORECASE)),
    3: (2, re.compile("U[1-4S]", re.IGNORECASE)),
    6: (5, re.compile(ident_chars, re.IGNORECASE)),
}
remapped_sorinds = {}
for sorind in sorind_list:
    if sorind == "Blank":
        continue

    matched = True
    matches_standard = re.match(naming_standard, sorind)
    subset_sorind, source_name = sorind[:12], sorind[12:]
    corrected = ''

    if matches_standard is None:
        matched = False
        for stind, descr in index_breakdown.items():
            enind, frmt = descr
            section = subset_sorind[stind:stind + enind]
            section_match = re.match(frmt, section)
            if section_match is None:
                if stind == 6 and section.lower() not in ident_breakdown.values():
                    # 'nsurf' and 'survy' should be 'graph'
                    corrected += f"{ident_breakdown['g']},"
            elif stind in (0, 3):
                section = section.upper()
                corrected += f"{section},"
            elif stind == 6:
                letter = section[0].lower()
                corrected += f"{ident_breakdown[letter]},"
        if corrected != 'Blank':
            corrected += source_name
    # There is a duplicate sorind in US4DE11M for Chart 12214
    if sorind == 'US,U1,graph,Chart 12214':
        corrected = 'US,US,graph,Chart 12214'
        matched = False
    # duplicate "graph" or "reprt" in the name -- found a graph,graph and graph,reprt
    # If not handled Bruty will report "Two IDs reference the same file ID"
    if re.search(r",\s*(graph|reprt),\s*(graph|reprt)", sorind, re.IGNORECASE):
        if corrected:
            pass  # I think this corrects the leading part of the sorind
        else:
            # remove all ",graph" or ",reprt" that are followed another ",graph" or ",reprt"
            # FIXME We have seen cases where the second of the two graph,reprt is correct so remove the first
            #   if the first occurence is correct we need better logic to find the matches
            corrected = re.sub(r"(,\s*(graph|reprt))(?=,\s*(graph|reprt))", "", sorind, flags=re.IGNORECASE)
            matched = False
    # double comma where only one should be -- US5ILMEC / DD-38023 had this
    # FIXME if sorind has two errors then only one will get corrected currently
    #   (i.e. corrected is not checked but is replaced)
    if re.search(r"\s*(graph|reprt),\s*,", sorind, re.IGNORECASE):
        if corrected:
            pass  # I think this corrects the leading part of the sorind
        else:
            # remove all ",graph" or ",reprt" that are followed another ",graph" or ",reprt"
            # FIXME We have seen cases where the second of the two graph,reprt is correct so remove the first
            #   if the first occurence is correct we need better logic to find the matches
            corrected = re.sub(r",\s*,", ",", sorind, flags=re.IGNORECASE)
            matched = False
    split_sorind = sorind.split(',')
    if bool(re.search(r'^\s+|\s+$', split_sorind[-1])):
        matched = False
        if corrected != '':
            split_sorind = corrected.split(',')
            corrected = (',').join([section.strip() for section in split_sorind])
        else:
            corrected = (',').join([section.strip() for section in split_sorind])

    if 'chart' in source_name:
        matched = False
        if corrected != '':
            corrected = corrected.replace('chart', 'Chart')
        else:
            corrected = sorind.replace('chart', 'Chart')

    if not matched:
        if corrected == 'Blank':
            LOGGER.warning(f"SORIND '{sorind}' could not be resolved to naming standards, replacing with: '{corrected}'")
        else:
            LOGGER.warning(f"SORIND '{sorind}' does not match naming standards, replacing with: '{corrected}'")
        if corrected not in remapped_sorinds:
            remapped_sorinds[corrected] = [sorind]
        else:
            remapped_sorinds[corrected].append(sorind)
    print(sorind, "   mapped to   "+corrected if corrected else " is ok")