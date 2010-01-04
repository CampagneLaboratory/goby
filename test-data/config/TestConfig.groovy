settings {
    // Values that are the same for ALL environments
    one = 1
}

environments {
    dubble {
        settings {
            two = 2 * 2
            three = "six"
            four = 8
            any = { String value -> "two-" + value }
        }
    }
    quadruple_four {
        settings {
            // For testing the merge
            four = "sixteen"
        }
    }
}
