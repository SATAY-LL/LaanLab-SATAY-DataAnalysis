def chromosomename_roman_to_arabic():
    """
    Creates two dictionaries for translating the chromosome names from roman to arabic numerals or vice versa.
    
    Syntax: arabic_to_roman, roman_to_arabic = chromosomename_roman_to_arabic()

    Returns
    -------

    arabic_to_roman : dict
    roman_to_arabic : dict

    """
    num_arabic = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]
    num_roman = [
        "I",
        "II",
        "III",
        "IV",
        "V",
        "VI",
        "VII",
        "VIII",
        "IX",
        "X",
        "XI",
        "XII",
        "XIII",
        "XIV",
        "XV",
        "XVI",
        "Mito",
    ]

    arabic_to_roman = dict(zip(num_arabic, num_roman))
    roman_to_arabic = dict(zip(num_roman, num_arabic))

    return arabic_to_roman, roman_to_arabic
