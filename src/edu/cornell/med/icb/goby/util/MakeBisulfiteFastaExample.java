/*
 * Copyright (C) 2009-2011 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.cornell.med.icb.goby.util;

import it.unimi.dsi.lang.MutableString;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

/**
 * @author Fabien Campagne
 *         Date: 11/25/11
 *         Time: 10:11 AM
 */
public class MakeBisulfiteFastaExample {
    static String[] example1_forward = {
            "m|995-995999-9399-9-99-9999999-99993-                       4                                 ",
            "r|GCCACCGGGCTGCCGAGTGGTCCCCGCCACCCCCACAGGCTCAGAAATTGTTCTCTCTGAAACTCTGAGCACCACCGCCCAGCCGGGAGAGG",
            "+|GTTATTGGGTTGCTGAGTGGTTTTTGTTATTTTTATAGGTTTAGAAATTGTTTTTTTTGAAATTTTGAGTATTATTGTTTAGTTGGGAGAGG",
            "+|GTTATTGGGTTGTTGAGTGGTTTTTGTTATTTTTATAGGTTTAGAAATTGTTTTTTTTGAAATTTTGAGTATTATTGTTTAGTTGGGAGAGG",
            "+|GTTATTGGGTTGTCGAGTGGTTTTTGTTATTTTTATAGGTTTAGAAATTGTTTTTTTTGAAATTTTGAGTATTATTGTTTAGTTGGGAGAGG",
            "+|GTTATTGGGTTGTCGAGTGGTTTTTGTTATTTTCATAGGTTTAGAAATTGTTTTTTTTGAAATTTTGAGTATTATTGTTTAGTTGGGAGAGG",
            "+|GTTATTGGGTTGTTGAGTGGTTTTTGTTATTTTCATAGGTTTAGAAATTGTTTTTTTTGAAATTTTGAGTATTATTGTTTAGTTGGGAGAGG",
            "+|GTCATTGGGTTGTTGAGTGGTTTTTGTTATTTTCATAGGTTTAGAAATTGTTTTTTTTGAAATTTTGAGTATTATTGTTTAGTTGGGAGAGG",
            "+|GTCATTGGGTTGTCGAGTGGTTTTTGTTATTTTTATAGGTTTAGAAATTGTTTTTTTTGAAATTTTGAGTATTATTGTTTAGTTGGGAGAGG",
            "+|GTCATTGGGTTGTTGAGTGGTTTTTGTTATTTTTATAGGTTTAGAAATTGTTTTTTTTGAAATTTTGAGTATTATTGTTTAGTTGGGAGAGG",
            "+|GTCATTGGGTTGTTGAGTGGTTTTTGTTATTTTTATAGGTTTAGAAATTGTTTTTTTTGAAATTTTGAGTATTATTGTTTAGTTGGGAGAGG",
            "+|GTCATTGGGTTGTTGAGTGGTTTTTGTTATTTTTATAGGTTTAGAAATTGTTTTTTTTGAAATTTTGAGTATTATTGTTTAGTTGGGAGAGG"};
    static String[] example1_reverse = {
            "m|               8                                3                            5              ",
            "r|CCTCTCCCGGCTGGGCGGTGGTGCTCAGAGTTTCAGAGAGAACAATTTCTGAGCCTGTGGGGGTGGCGGGGACCACTCGGCAGCCCGGTGGC",
            "-|TTTTTTTTGGTTGGGCGGTGGTGTTTAGAGTTTTAGAGAGAATAATTTCTGAGTTTGTGGGGGTGGTGGGGATTATTTGGTAGTTTGGTGGT",
            "-|TTTTTTTTGGTTGGGCGGTGGTGTTTAGAGTTTTAGAGAGAATAATTTCTGAGTTTGTGGGGGTGGTGGGGATTATTTGGTAGTTTGGTGGT",
            "-|TTTTTTTTGGTTGGGCGGTGGTGTTTAGAGTTTTAGAGAGAATAATTTCTGAGTTTGTGGGGGTGGTGGGGATTATTCGGTAGTTTGGTGGT",
            "-|TTTTTTTTGGTTGGGCGGTGGTGTTTAGAGTTTTAGAGAGAATAATTTTTGAGTTTGTGGGGGTGGTGGGGATTATTCGGTAGTTTGGTGGT",
            "-|TTTTTTTTGGTTGGGCGGTGGTGTTTAGAGTTTTAGAGAGAATAATTTTTGAGTTTGTGGGGGTGGTGGGGATTATTCGGTAGTTTGGTGGT",
            "-|TTTTTTTTGGTTGGGCGGTGGTGTTTAGAGTTTTAGAGAGAATAATTTTTGAGTTTGTGGGGGTGGTGGGGATTATTCGGTAGTTTGGTGGT",
            "-|TTTTTTTTGGTTGGGCGGTGGTGTTTAGAGTTTTAGAGAGAATAATTTTTGAGTTTGTGGGGGTGGTGGGGATTATTCGGTAGTTTGGTGGT",
            "-|TTTTTTTTGGTTGGGCGGTGGTGTTTAGAGTTTTAGAGAGAATAATTTTTGAGTTTGTGGGGGTGGTGGGGATTATTTGGTAGTTTGGTGGT",
            "-|TTTTTTTTGGTTGGGTGGTGGTGTTTAGAGTTTTAGAGAGAATAATTTTTGAGTTTGTGGGGGTGGTGGGGATTATTTGGTAGTTTGGTGGT",
            "-|TTTTTTTTGGTTGGGTGGTGGTGTTTAGAGTTTTAGAGAGAATAATTTTTGAGTTTGTGGGGGTGGTGGGGATTATTTGGTAGTTTGGTGGT"
    } ;

    public static void main(String args[]) throws IOException {

        PrintWriter writer = new
                PrintWriter(System.out);

        convert(writer, example1_forward);
        convert(writer, example1_reverse);
        writer.flush();
    }

    private static void convert(PrintWriter writer, String[] example) {
        ArrayList<String> reads = new ArrayList(filterReads(example));
        int index = 1;
        for (String r : reads) {
            writer.printf(">%d%c %n%s%n", index++, r.charAt(0), r.substring(2));
        }
    }

    private static ArrayList<String> filterReads(String[] example) {
        ArrayList<String> result = new ArrayList<String>();
        for (String s : example) {
            if (s.startsWith("+|")) {
                result.add(s);
            } else {
                if (s.startsWith("-|")) {
                    /*// reverse read sequence
                    MutableString reverse = new MutableString();
                    reverse.append("-|");
                    for (int i = s.length() - 1; i >= 2; i--) {
                        reverse.append(s.charAt(i));
                    }
                    result.add(reverse.toString());      */
                    result.add(s);
                }
            }

        }
        return result;
    }


    static String example2 = "                 reference (null)/0 GCCACCGGGCTGCCGAGTGGTCCCCGCCACCCCCACAGGCTCAGAAATTGTTCTCTCTGAAACTCTGAGCACCACCGCCCAGCCGGGAGAGGCCAGGGTCTGCAGCCCAGCGCCTGAAGC\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/109866693 ---....T..T........TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/140194384 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/62615418 ---T...T..T........TTT..TT.TTTTT.T...T.T..........T.T.\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/139576185 ---T...T..T........TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/36635724 ---T...T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/130433408 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/95254095 ---....T..T........TTT..TT.TTTTT.T...T.T..........T.T.\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/104422232 ---T...T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/94619988 ---....T..T........TTT..TT.TTTTT.T...T.T..........T.T.\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/105051398 ---....T..T........TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/95088395 ---T...T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/84563696 ---T...T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/79659300 ---T...T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/85148239 ---....T..T........TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/77375093 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/67391393 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/76981077 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/67467764 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/89913671 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/25048859 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/44962489 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/132738042 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/46766308 ---....T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/133186407 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/46863224 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/126875344 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/47357085 ---....T..T........TTT..TT.TTTTT.T...T.T..........T.T.\n" +
            "         KIJECNH-DLBCL-nov-21-NB500-NB-_-_/9547675 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/16733737 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "         KIJECNH-DLBCL-nov-21-NB500-NB-_-_/3877656 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/14693947 ---....T..T........TTT..TT.TTTTT.T...T.T..........T.T.\n" +
            "         KIJECNH-DLBCL-nov-21-NB500-NB-_-_/4668429 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/34980515 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "         KIJECNH-DLBCL-nov-21-NB500-NB-_-_/4421260 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/34635132 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/31781300 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/102204701 ---....T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/31872226 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/102224110 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/32386218 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/102373600 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/59760433 ---T...T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/49147994 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/91985010 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/48986218 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/26580976 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/75221014 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/54189220 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/57080849 ---T...T..T........TTT..TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/49396990 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/41885477 ---T...T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/131193858 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/73889240 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/43262875 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/73206822 ---....T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/43363445 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/78722440 ---T...T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/42788880 ---....T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/77907932 ---T...T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/85417259 ---T...T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/20976688 ---T...T..TT.......TTTT.TT.TTTTT.TG..T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/75934353 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/55431080 ---....T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/76252971 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/55279591 ---T...T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/105944383 ---T...T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/55686813 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/105821843 ---....T..T........TTT..TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/56158541 ---....T..TT.......TTT..TTGTTTTT.T...T.T..........T.T.\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/118534480 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/89116352 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/136903805 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/88036447 ---T...T..T........TTT..TT.TTTTT.T...T.T..........T.T.\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/136014994 ---T...T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/93294429 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/136049357 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/93415057 ---T...T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/53174831 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/93009372 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/50159697 ---T...T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/10384638 ---T...T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/50339528 ---T...T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/10399836 ---....T..T........TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/50808302 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/36221913 ---T...T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.\n" +
            "          KIJECNH-DLBCL-nov-21-NB500-NB-_-_/202164 ---T...T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/126219326 ---T...T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.\n" +
            "          KIJECNH-DLBCL-nov-21-NB500-NB-_-_/338184 ---T...T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/126291148 ---T...T..T........TTT..TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/83423603 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/125723040 ---T...T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/82860476 ---....T..T........TTT..TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/48422155 ---T...T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/121571750 ---T...T..TT.......TTTT.TT.GTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/48116969 ---....T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/101466391 ---T...T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/80788831 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/31327496 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/63760895 ---T...T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/103990042 ---....T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/58820462 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/80296654 ---T...T..T........TTT..TT.TTTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/68193725 ---T...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.\n" +
            "         KIJECNH-DLBCL-nov-21-NB500-NB-_-_/8161605 ----...T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.T\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/130823054 ----...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.T\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/53355143 ----...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.T\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/128699950 ----...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.T\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/115412230 ----...T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.T\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/108540496 ----...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.T\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/28291940 ----...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.T\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/13644018 ----...T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.T\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/83612791 ----...T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.T\n" +
            "          KIJECNH-DLBCL-nov-21-NB500-NB-_-_/967891 ----...T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.T\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/60798800 ----...T..T........TTTT.TT.TTTTT.T...T.T..........T.T.T\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/131345168 ----...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.T\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/115418146 ----...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.T\n" +
            "          KIJECNH-DLBCL-nov-21-NB500-NB-_-_/962045 ----...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.T\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/36978546 ----...T..T........TTT..TT.TTTTT.T...T.T..........T.T.T\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/94961389 ----...T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.T\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/107617902 ----...T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.T\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/77039706 ----...T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.T\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/45160172 ----...T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.T\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/69426227 ----...T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.T\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/24068828 ----...T..TT.......TTTT.TT.TTTTT.T...T.T..........T.T.T\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/15054557 ----...T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.T\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/31866669 ----...T..TT.......TTT..TT.TTTTT.T...T.T..........T.T.T\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/116863923 -------T..TT.......TTT..TT.ATTTT.T...T.T..........T.T.\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/22137777 --------..TT.......TTTT.TT.TGTTT.T...T.T..........T.T.T\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/116801860 --------..TT.......TTT..TT.TGTTT.T...T.T..........T.T.T\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/139683787 -------------------------------....AA....A.....A........A.......A.A.......A....A..\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/61927158 -------------------------------....AA....A.....A........A.......A.A.......A....A..\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/49998120 -------------------------------....AA....A.....A........A.......A.A............A..\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/51593163 -------------------------------....AA....A.....A........A.......A.A............A..\n" +
            "         KIJECNH-DLBCL-nov-21-NB500-NB-_-_/9536119 -------------------------------....AA....A.....A........A.......A.A.......A....A..\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/59088148 -------------------------------....AA....A.....A........A.......A.A.......A....A..\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/74213699 -------------------------------....AA....A.....A........A.......A.A.......A....A..\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/53800289 -------------------------------....AA....A.....A........A.......A.A.......A....A..\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/51968806 -------------------------------....AA....A.....A........A.......A.A.......A....A..\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/75907100 -------------------------------....AA....A.....A........A.......A.A.......A....A..\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/91469889 -------------------------------....AA....A.....A........A.......A.A.......A....A..\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/65613830 -------------------------------....AA....A.....A........A.......A.A.......A....A..\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/20981877 -------------------------------....AA....A.....A........A.......A.A.......A....A..\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/85975009 -------------------------------....AA....A.....A........A.......A.A.......A....A..\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/118004751 -------------------------------....AA....A.....A........A.......A.A.......A....A..\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/55169664 -------------------------------....AA....A.....A........A.......A.A.......A....A..\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/106401017 -------------------------------....AA....A.....A........A.......A.A.......A....A..\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/33486465 -------------------------------....AA....A.....A........A.......A.A.......A....A..\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/134153282 -------------------------------....AA....A.....A........A.......A.A.......A....A..\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/125660743 -------------------------------....AA....A.....A........A.......A.A.......A....A..\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/134420894 --------------------------------...AA....A.....A........A.......A.A............A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/56278265 --------------------------------...AA....A.....A.C......A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/90560120 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/41340703 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/90886945 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/41243361 --------------------------------...AA....A.....A........A.......A.A............A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/91263804 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/16114371 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/91099908 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/113414497 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/91484282 --------------------------------...AA....A.....A........A......CA.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/81354840 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/124219359 --------------------------------...AA....A.....A........A.......A.........A....A...\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/80789455 --------------------------------...AA....A.....A........A.......A.A.......A....A...\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/124144693 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/81554832 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/123746836 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/75690596 --------------------------------...AA....A.....A........A.......A.A............A...\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/123670810 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/75732355 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/93764389 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/75766858 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/45949956 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/43804723 --------------------------------...AA....A.....A.C......A.C.....A.A.......A....\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/78354517 --------------------------------...AA....A.....A........A.......A.A............A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/43804724 --------------------------------...AA....A.....A.C......A.CG....A.A.......A....\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/118953967 --------------------------------...AA....A.....AC.......A....C..A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/96681992 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/118984905 --------------------------------...AA....A.....A........A....A..A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/96365446 --------------------------------...AA....A.....A........A.......A.A............A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/118239264 --------------------------------...AA....A.....A........A.......A.A............A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/96417731 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "         KIJECNH-DLBCL-nov-21-NB500-NB-_-_/5562681 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/65773713 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/25862825 --------------------------------...AA....A.....A........A.......A.A............A...\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/66317535 --------------------------------...AA....A.....A........A.......A.A.......A....A...\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/104081410 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/68451557 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/103264308 --------------------------------...AA....A.....A........A.......A.A............A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/136368166 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/23170736 --------------------------------...AA....A.....A........A.......A.A.......A....\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/136925977 --------------------------------...AA....A.....A........A.......A.A.......A....A...\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/88156259 --------------------------------...AA....A.....A........A.......A.........A....A...\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/70986713 --------------------------------...AA....A.....A........A.......A.A.C..........A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/106036495 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/71439880 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/106300685 --------------------------------...AA....A.....A........A.......A.A............A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/27849693 --------------------------------...AA....A.....A........A.......A.A............A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/106701183 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/28546609 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/58725558 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/10517153 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/57778511 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/48283630 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/57823780 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/128524304 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/57619316 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/128974009 --------------------------------...AA....A.....A........A.......A.A............A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/98147420 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/36346940 --------------------------------...AA....A.....A........A.......A.A............A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/98498664 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "         KIJECNH-DLBCL-nov-21-NB500-NB-_-_/3449998 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/51302168 --------------------------------...AA....A.....A........A.......A.A.......A....A...\n" +
            "         KIJECNH-DLBCL-nov-21-NB500-NB-_-_/3467733 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/50941679 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/120707624 --------------------------------...AA....A.....A........A.C.....A.A............\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/50795960 --------------------------------...AA....A.....A........A.......A.A............A...\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/121362740 --------------------------------...AA....A.....A........A.C.....A.A.......A....\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/31363926 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/13583660 --------------------------------...AA....A.....A........A.......A.A............A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/130915255 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/115578984 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/121258350 --------------------------------...AA....A.....A........A.C.....A.A.......A....\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/30570628 --------------------------------...AA....C.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/13175092 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/116153756 --------------------------------...AA....A.....A........A.......A.A............A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/38371115 --------------------------------...AA....A.....A........A.......A.A............A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/116174405 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/62606042 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/31772036 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/110086527 --------------------------------...AA....A.....A........A.......A.A.......A....A...\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/31981244 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/39590555 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/19392704 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/40049027 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/19081182 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/39703481 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/45073856 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/39721483 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/44849691 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/11653466 --------------------------------...AA....A.....A........A.......A.A............A...\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/79733697 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/11739943 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/22209968 --------------------------------...AA....A.....A........A.C.....A.A.......A.......A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/12015941 --------------------------------...AA....A.....A........A.......A.A............A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/21434189 --------------------------------...AA....A.....A........A.......A.A............A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/127694003 --------------------------------...AA....A.....A........A.......A.A............A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/130180350 --------------------------------...AA....A.....A........A.......A.A............A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/126921670 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/129859814 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/81658572 --------------------------------...AA....A.....A........A.......A.A............A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/129748770 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/82310016 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/105453651 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/82702420 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/104689690 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/137956266 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/105164401 --------------------------------...AA....A.....AC.......A.......A.A.......A....A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/137052408 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/105228468 --------------------------------...AA....A.....A........A.......A.A............A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/122861641 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/104956351 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/121725192 --------------------------------...AA....A.....A........A.......A.A............\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/99393418 --------------------------------...AA....A.....A........A.......A.A.......A....A...\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/121756958 --------------------------------...AA....A.....A........A.T.....A.A.......A....\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/89872822 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/95178385 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/89959198 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/94673634 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/125081379 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/112417528 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/125242886 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/25028400 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/14435252 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/47220593 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/64304309 --------------------------------...AA....A.....A........A.......A.A............A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/46875109 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/65130594 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/46891988 --------------------------------...AA....A.....A........A.......A.A............A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/64046244 --------------------------------...AA....A.....A........A.......A.A............A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/47121289 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/64078522 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/69464652 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "         KIJECNH-DLBCL-nov-21-NB500-NB-_-_/4959998 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/17099849 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "         KIJECNH-DLBCL-nov-21-NB500-NB-_-_/4054484 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/16828421 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/56946008 --------------------------------...AA....A.....A........A.......A.A............A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/35142837 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/120259113 --------------------------------...AA....A.....A........A.C.....A.A.......A....\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/34892938 --------------------------------...AA....A.....A........A.......A.A............A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/120029517 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/34286579 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/135321554 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/102032545 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "         KIJECNH-DLBCL-nov-21-NB500-NB-_-_/6592668 --------------------------------...AA....A.....A........A.......A.A............A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/102880594 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/114676928 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/58985078 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "         KIJECNH-DLBCL-nov-21-NB500-NB-_-_/9861923 --------------------------------...AA....A.....A........A.......A.A.......A....A...\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/49948574 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "       KIJECNH-DLBCL-nov-21-NB500-NB-_-_/114712605 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/26812857 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "         KIJECNH-DLBCL-nov-21-NB500-NB-_-_/9129764 --------------------------------...AA....A.....A.......CA.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/72316399 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/91933466 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/27141839 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "         KIJECNH-DLBCL-nov-21-NB500-NB-_-_/9553609 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/27399914 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "         KIJECNH-DLBCL-nov-21-NB500-NB-_-_/9681328 --------------------------------...AA....A.....A........A.......A.A.......A....A..A\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/77764790 ---------------------------------..AA....A.....A........A.......A.A.......A....\n" +
            "        KIJECNH-DLBCL-nov-21-NB500-NB-_-_/40978275 ----------------------------------------.A.....A........A.......A.A.......A....A..A";
}
