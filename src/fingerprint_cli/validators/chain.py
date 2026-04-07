from .common import check_exists, check_tool

def validate_chain_inputs(old_ref, new_ref):
    check_exists(old_ref, "Old reference genome")
    check_exists(new_ref, "New reference genome")

    for tool in ["minimap2", "transanno"]:
        check_tool(tool)
