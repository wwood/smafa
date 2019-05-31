

fn translate(x: &[u8; 3]) -> u8 {
    match x {
        [b'T',b'T',b'T'] =>	b'F',
        [b'T',b'T',b'C'] =>	b'F',
        [b'T',b'T',b'A'] =>	b'L',
        [b'T',b'T',b'G'] =>	b'L',
        [b'C',b'T',b'T'] =>	b'L',
        [b'C',b'T',b'C'] =>	b'L',
        [b'C',b'T',b'A'] =>	b'L',
        [b'C',b'T',b'G'] =>	b'L',
        [b'A',b'T',b'T'] =>	b'I',
        [b'A',b'T',b'C'] =>	b'I',
        [b'A',b'T',b'A'] =>	b'I',
        [b'A',b'T',b'G'] =>	b'M',
        [b'G',b'T',b'T'] =>	b'V',
        [b'G',b'T',b'C'] =>	b'V',
        [b'G',b'T',b'A'] =>	b'V',
        [b'G',b'T',b'G'] =>	b'V',
        [b'T',b'C',b'T'] =>	b'S',
        [b'T',b'C',b'C'] =>	b'S',
        [b'T',b'C',b'A'] =>	b'S',
        [b'T',b'C',b'G'] =>	b'S',
        [b'C',b'C',b'T'] =>	b'P',
        [b'C',b'C',b'C'] =>	b'P',
        [b'C',b'C',b'A'] =>	b'P',
        [b'C',b'C',b'G'] =>	b'P',
        [b'A',b'C',b'T'] =>	b'T',
        [b'A',b'C',b'C'] =>	b'T',
        [b'A',b'C',b'A'] =>	b'T',
        [b'A',b'C',b'G'] =>	b'T',
        [b'G',b'C',b'T'] =>	b'A',
        [b'G',b'C',b'C'] =>	b'A',
        [b'G',b'C',b'A'] =>	b'A',
        [b'G',b'C',b'G'] =>	b'A',
        [b'T',b'A',b'T'] =>	b'Y',
        [b'T',b'A',b'C'] =>	b'Y',
        [b'T',b'A',b'A'] =>	b'*',
        [b'T',b'A',b'G'] =>	b'*',
        [b'C',b'A',b'T'] =>	b'H',
        [b'C',b'A',b'C'] =>	b'H',
        [b'C',b'A',b'A'] =>	b'Q',
        [b'C',b'A',b'G'] =>	b'Q',
        [b'A',b'A',b'T'] =>	b'N',
        [b'A',b'A',b'C'] =>	b'N',
        [b'A',b'A',b'A'] =>	b'K',
        [b'A',b'A',b'G'] =>	b'K',
        [b'G',b'A',b'T'] =>	b'D',
        [b'G',b'A',b'C'] =>	b'D',
        [b'G',b'A',b'A'] =>	b'E',
        [b'G',b'A',b'G'] =>	b'E',
        [b'T',b'G',b'T'] =>	b'C',
        [b'T',b'G',b'C'] =>	b'C',
        [b'T',b'G',b'A'] =>	b'*',
        [b'T',b'G',b'G'] =>	b'W',
        [b'C',b'G',b'T'] =>	b'R',
        [b'C',b'G',b'C'] =>	b'R',
        [b'C',b'G',b'A'] =>	b'R',
        [b'C',b'G',b'G'] =>	b'R',
        [b'A',b'G',b'T'] =>	b'S',
        [b'A',b'G',b'C'] =>	b'S',
        [b'A',b'G',b'A'] =>	b'R',
        [b'A',b'G',b'G'] =>	b'R',
        [b'G',b'G',b'T'] =>	b'G',
        [b'G',b'G',b'C'] =>	b'G',
        [b'G',b'G',b'A'] =>	b'G',
        [b'G',b'G',b'G'] =>	b'G',
        _ => 'X' as u8,
    }
}

pub fn translate_codons(input_dna_sequence: &[u8]) -> Vec<u8> {
    let mut to_return:Vec<u8> = vec!();
    let mut i: usize = 0;
    while i < input_dna_sequence.len() {
        to_return.push(translate(&[
            input_dna_sequence[i],
            input_dna_sequence[i+1],
            input_dna_sequence[i+2]
        ]));
        i += 3;
        if i > input_dna_sequence.len() {
            panic!("Sequence length not divisible by 3");
        }
    }
    return to_return;
}


#[cfg(test)]
mod tests {
    use translation::translate_codons;

    #[test]
    fn test_hello_world() {
        assert_eq!(vec!(b'F',b'L'), translate_codons(b"TTTTTA"))
    }
}
