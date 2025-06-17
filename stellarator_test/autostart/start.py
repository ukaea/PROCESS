from stellarator_test.autostart import generate_input, run_cases, collect_results

case_name = 'low_blanket'
prefix = 'squid'

generate_input.main(case_name, prefix = prefix)
run_cases.main(case_name, prefix = prefix)
collect_results.main(case_name, prefix = prefix)
