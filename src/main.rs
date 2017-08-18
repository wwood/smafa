extern crate bio;

use std::str;

use bio::alphabets::dna;
use bio::data_structures::fmindex::{FMIndex};
use bio::data_structures::suffix_array::{suffix_array};
use bio::data_structures::bwt::{bwt, less, Occ};
use bio::data_structures::fmindex::FMIndexable;

use std::io;
use bio::io::fasta;

use std::time::Instant;

fn main(){
    let mut text = vec![];

    let reader = fasta::Reader::new(io::stdin());

    let mut now = Instant::now();
    let mut new_now;

    {
        let mut i: u64 = 0;
        for record in reader.records() {
            text.extend(record.unwrap().seq().iter().cloned());
            text.extend_from_slice(b"$");
            i+=1;
            if i % 10000 == 0 {
                print!(".")
            }
        }
        println!();
    }

    // println!("total: {:?}", text)
    let alphabet = dna::n_alphabet();
    new_now = Instant::now(); println!("reading {:?}", new_now.duration_since(now)); now = new_now;
    let sa = suffix_array(&text);
    new_now = Instant::now(); println!("suffix array {:?}", new_now.duration_since(now)); now = new_now;
    let bwt = bwt(&text, &sa);
    new_now = Instant::now(); println!("bwt {:?}", new_now.duration_since(now)); now = new_now;
    let less = less(&bwt, &alphabet);
    new_now = Instant::now(); println!("less {:?}", new_now.duration_since(now)); now = new_now;
    let occ = Occ::new(&bwt, 3, &alphabet);
    new_now = Instant::now(); println!("occ {:?}", new_now.duration_since(now)); now = new_now;
    let fm = FMIndex::new(&bwt, &less, &occ);
    new_now = Instant::now(); println!("fmindex {:?}", new_now.duration_since(now)); now = new_now;

    let pattern = b"ATGC";
    //let intervals = fmdindex.smems(pattern, 2);
    let intervals = fm.backward_search(pattern.iter());
    new_now = Instant::now(); println!("search {:?}", new_now.duration_since(now)); now = new_now;
    println!("occ2..");
    // println!("fwd: {:?}", intervals);
    let some = intervals.occ(&sa);
    new_now = Instant::now(); println!("occ2 {:?}", new_now.duration_since(now));// now = new_now;
    println!("found: {:} hits", some.len());
}

// fn main2(){

//     // let text = b"ATCCC$GGGAT$AGGGA$";
//     // let alphabet = dna::n_alphabet();
//     // let sa = suffix_array(text);
//     // let bwt = bwt(text, &sa);
//     // let less = less(&bwt, &alphabet);
//     // let occ = Occ::new(&bwt, 3, &alphabet);
//     // let fm = FMIndex::new(&bwt, &less, &occ);
//     // //let fmdindex = FMDIndex::from(fm);

//     // let pattern = b"GGGA";
//     // //let intervals = fmdindex.smems(pattern, 2);
//     // let intervals = fm.backward_search(pattern.iter());
//     // println!("fwd: {:?}", intervals);
//     // let positions = intervals.occ(&sa);
//     // println!("fwd: {:?}", positions);
//     // // let forward_positions = intervals[0].forward().occ(&sa);

//     // println!("fwd: {:?}", forward_positions);


//     // let text = b"ATCCC$GGGAT$AGGGA$";
//     // let alphabet = dna::n_alphabet();
//     // let sa = suffix_array(text);
//     // let bwt = bwt(text, &sa);
//     // let less = less(&bwt, &alphabet);
//     // let occ = Occ::new(&bwt, 3, &alphabet);
//     // let fm = FMIndex::new(&bwt, &less, &occ);
//     // search(fm, sa, b"GG")
// }

// pub fn search(fm_index: &FMIndex<DerefBWT<Target=BWT>, DerefLess<Target=Less>, DerefOcc<Target=Occ>>, sa: &SuffixArray, search_pattern: [u8]){
//     let intervals = fm_index.backward_search(search_pattern.iter());
//     println!("fwd: {:?}", intervals.occ(sa));
// }
