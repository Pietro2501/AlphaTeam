def parse_fastq(file_obj, offset=33) -> dict:
    """
    Parses a FASTQ file and returns a dictionary representation of the sequences,
    quality scores, and additional metadata. The function processes the input file
    line by line, ensuring the structural integrity for each FASTQ record, including
    checking coherence of headers and applicable offset values for quality score
    conversion.

    :param file_obj: A file object containing the FASTQ data to be parsed.
    :type file_obj: str
    :param offset: Offset value used for ASCII quality score conversion.
                   Should only be either 33 or 64.
    :type offset: int
    :return: A dictionary where the keys are sequence identifiers (headers starting
             with '@') and values are sub-dictionaries containing sequence strings,
             ASCII quality scores, and numeric quality scores. Returns None if the file
             is found to be corrupted or improperly formatted.
    :rtype: dict or None
    """
    with open(file_obj, "r") as file_obj:
        def check_offset(offset: int) -> bool:
            res = False
            if type(offset) == int and offset in [33, 64]:
                res = True
            return res

        def check_coherence(acc: str, plus: str) -> bool:
            """
            Validates the coherence of a FASTQ record by checking the format of
            the accession line and the "+" line.

            The function ensures that the accession line starts with "@" and
            the "+" line starts with "+", which is required for the correct
            format of a FASTQ record.

            :param acc: The accession line from the FASTQ record.
            :type acc: str
            :param plus: The "+" line from the FASTQ record.
            :type plus: str
            :return: True if both the accession and the "+" lines are in the
                correct format, otherwise False.
            :rtype: bool
            """
            res = False
            if acc.startswith("@") and plus.startswith("+"):
                res = True
            return res

        if check_offset(offset):
            fastq_dict = {}
            line = file_obj.readline()
            while line:
                acc = line.strip()
                fastq_dict.setdefault(acc, {})
                # qui andrebbe inserito un try/expcet per verificare che effettivamente tutte le linee siano lette
                seq, plus, qual = file_obj.readline().strip(), file_obj.readline().strip(), file_obj.readline().strip()
                if not check_coherence(acc, plus):
                    print("Il file è corrotto!!!")
                    fastq_dict = None
                    break
                if seq and qual:
                    fastq_dict[acc]["seq"] = seq
                    fastq_dict[acc]["ASCII_qual"] = qual
                    fastq_dict[acc]["qual"] = [ord(q) - offset for q in fastq_dict[acc]["ASCII_qual"]]
                    # fastq_dict[acc]["Pe"] = [10 ** (q / -10) for q in fastq_dict[acc]["qual"]]
                else:
                    print("Il file è corrotto!!!")
                    fastq_dict = None
                    break
                line = file_obj.readline()
            return fastq_dict
        elif not check_offset(offset):
            print(f"l'offset indicato non è numerico o è un valore differente da 33 o 64!!!")