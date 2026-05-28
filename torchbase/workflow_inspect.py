"""Workflow inspection and visualization for WDL workflows."""

import re
from pathlib import Path
from typing import Optional, Dict, Tuple


class WDLParser:
    """Simple WDL parser for extracting task flow and structure."""

    def __init__(self, wdl_content: str):
        """Initialize parser with WDL content."""
        self.content = wdl_content
        self.tasks = {}
        self.workflow_name = None
        self.workflow_inputs = {}
        self.workflow_outputs = {}
        self.task_calls = []
        self.conditionals = []
        self.parse()

    def parse(self):
        """Parse WDL content to extract workflow structure."""
        # First, validate basic WDL syntax
        self._validate_syntax()

        try:
            # Extract workflow name
            workflow_match = re.search(r'workflow\s+(\w+)', self.content)
            if workflow_match:
                self.workflow_name = workflow_match.group(1)

            # Extract workflow inputs
            self._extract_workflow_inputs()

            # Extract workflow outputs
            self._extract_workflow_outputs()

            # Extract task calls
            self._extract_task_calls()

            # Extract task definitions
            self._extract_tasks()

            # Extract conditionals
            self._extract_conditionals()
        except Exception as e:
            raise ValueError(f"Failed to parse WDL: {e}")

    def _validate_syntax(self):
        """Check for basic WDL syntax errors."""
        # Check for version declaration
        if not re.search(r'version\s+\d+\.\d+', self.content):
            raise ValueError("WDL must include version declaration (e.g., 'version 1.0')")

        # Check for unmatched braces/brackets
        open_braces = self.content.count('{')
        close_braces = self.content.count('}')
        if open_braces != close_braces:
            raise ValueError(f"WDL syntax error: mismatched braces ({open_braces} open, {close_braces} close)")

        # Note: We don't validate import paths here since they're relative to the WDL file location
        # and may not be resolvable from the current working directory. The WDL engine will
        # handle import resolution during actual execution.

        # Check for clearly malformed syntax patterns
        if re.search(r'\{\{\{', self.content) or re.search(r'\}\}\}', self.content):
            raise ValueError("WDL syntax error: malformed brace syntax")

    def _extract_workflow_inputs(self):
        """Extract workflow input section."""
        input_match = re.search(
            r'workflow\s+\w+\s*\{[^}]*?input\s*\{([^}]+)\}',
            self.content,
            re.DOTALL
        )
        if input_match:
            input_section = input_match.group(1)
            self.workflow_inputs = self._parse_declarations(input_section)

    def _extract_workflow_outputs(self):
        """Extract workflow output section."""
        output_match = re.search(
            r'output\s*\{([^}]+)\}',
            self.content,
            re.DOTALL
        )
        if output_match:
            output_section = output_match.group(1)
            for line in output_section.split('\n'):
                line = line.strip()
                if line and not line.startswith('//'):
                    # Parse output declaration
                    parts = re.match(r'(\w+)\s+(\w+)\s*=', line)
                    if parts:
                        type_name = parts.group(1)
                        var_name = parts.group(2)
                        self.workflow_outputs[var_name] = type_name

    def _extract_task_calls(self):
        """Extract call statements from workflow body."""
        # Find workflow body
        workflow_match = re.search(
            r'workflow\s+\w+\s*\{(.*?)\n\s*output\s*\{',
            self.content,
            re.DOTALL
        )
        if workflow_match:
            workflow_body = workflow_match.group(1)

            # Find all call statements
            call_pattern = r'call\s+(\w+)(?:\s+as\s+(\w+))?'
            for match in re.finditer(call_pattern, workflow_body):
                task_name = match.group(1)
                alias = match.group(2) or task_name
                self.task_calls.append({'task': task_name, 'alias': alias})

    def _extract_tasks(self):
        """Extract task definitions."""
        task_pattern = r'task\s+(\w+)\s*\{([^}]+?)(?=task\s+\w+|$)'
        for match in re.finditer(task_pattern, self.content, re.DOTALL):
            task_name = match.group(1)
            task_body = match.group(2)

            inputs = self._extract_task_inputs(task_body)
            outputs = self._extract_task_outputs(task_body)

            self.tasks[task_name] = {
                'name': task_name,
                'inputs': inputs,
                'outputs': outputs
            }

    def _extract_task_inputs(self, task_body: str) -> Dict[str, str]:
        """Extract inputs from a task body."""
        input_match = re.search(r'input\s*\{([^}]+)\}', task_body, re.DOTALL)
        if input_match:
            input_section = input_match.group(1)
            return self._parse_declarations(input_section)
        return {}

    def _extract_task_outputs(self, task_body: str) -> Dict[str, str]:
        """Extract outputs from a task body."""
        output_match = re.search(r'output\s*\{([^}]+)\}', task_body, re.DOTALL)
        if output_match:
            output_section = output_match.group(1)
            outputs = {}
            for line in output_section.split('\n'):
                line = line.strip()
                if line and not line.startswith('//'):
                    parts = re.match(r'(\w+)\s+(\w+)\s*=', line)
                    if parts:
                        type_name = parts.group(1)
                        var_name = parts.group(2)
                        outputs[var_name] = type_name
            return outputs
        return {}

    def _parse_declarations(self, section: str) -> Dict[str, Tuple[str, Optional[str]]]:
        """Parse variable declarations from a section.

        Returns dict mapping var_name to (type_name, default_value_or_none)
        """
        declarations = {}
        for line in section.split('\n'):
            line = line.strip()
            if line and not line.startswith('//'):
                # Match: Type name [= default]
                match = re.match(r'(\w+)\s+(\w+)(?:\s*=\s*(.+))?', line)
                if match:
                    type_name = match.group(1)
                    var_name = match.group(2)
                    default_val = match.group(3).strip() if match.group(3) else None
                    declarations[var_name] = (type_name, default_val)
        return declarations

    def _extract_conditionals(self):
        """Extract if statements from workflow."""
        if_pattern = r'if\s*\(([^)]+)\)\s*\{'
        for match in re.finditer(if_pattern, self.content):
            condition = match.group(1).strip()
            self.conditionals.append(condition)


class WorkflowDiagramRenderer:
    """Renders WDL workflow as ASCII box diagram."""

    def __init__(self, parser: WDLParser, verbose: bool = False):
        """Initialize renderer."""
        self.parser = parser
        self.verbose = verbose
        self.lines = []

    def render(self) -> str:
        """Render the workflow diagram."""
        self.lines = []

        # Header
        self._add_line("┌" + "─" * 78 + "┐")
        self._add_line("│ Workflow: " + (self.parser.workflow_name or "unknown").ljust(67) + "│")
        self._add_line("├" + "─" * 78 + "┤")

        # Inputs section
        if self.parser.workflow_inputs:
            self._add_line("│ Inputs:".ljust(80) + "│")
            for var_name, type_info in self.parser.workflow_inputs.items():
                if isinstance(type_info, tuple):
                    type_name, default_val = type_info
                else:
                    type_name, default_val = type_info, None

                opt = "?" if "?" in type_name else ""
                # Mark parameters with defaults as optional
                if default_val and not opt:
                    opt = "?"

                if self.verbose and default_val:
                    line = f"│   • {var_name}: {type_name}{opt} = {default_val}"
                else:
                    line = f"│   • {var_name}: {type_name}{opt}"
                self._add_line(line.ljust(80) + "│")
            self._add_line("├" + "─" * 78 + "┤")

        # Task calls section
        if self.parser.task_calls:
            self._add_line("│ Task Flow:".ljust(80) + "│")

            for i, call in enumerate(self.parser.task_calls):
                task_name = call['task']
                task_info = self.parser.tasks.get(task_name, {})

                # Connector line before task (except first)
                if i > 0:
                    self._add_line("│    ↓".ljust(80) + "│")

                # Task box
                self._add_line(f"│  ┌─ {task_name}".ljust(80) + "│")

                # Show key inputs
                if task_info.get('inputs') and not self.verbose:
                    # Show only File inputs
                    file_inputs = {}
                    for k, v in task_info['inputs'].items():
                        type_str = v[0] if isinstance(v, tuple) else v
                        if 'File' in type_str or 'Int' in type_str or 'Boolean' in type_str:
                            file_inputs[k] = v
                    for var_name, type_info in list(file_inputs.items())[:3]:
                        type_str = type_info[0] if isinstance(type_info, tuple) else type_info
                        line = f"│  │  {var_name}: {type_str}"
                        self._add_line(line.ljust(80) + "│")
                elif self.verbose and task_info.get('inputs'):
                    for var_name, type_info in task_info['inputs'].items():
                        if isinstance(type_info, tuple):
                            type_str, default_val = type_info
                            opt_indicator = "?" if "?" in type_str else ""
                            if default_val:
                                line = f"│  │  {var_name}: {type_str}{opt_indicator} = {default_val}"
                            else:
                                line = f"│  │  {var_name}: {type_str}{opt_indicator}"
                        else:
                            line = f"│  │  {var_name}: {type_info}"
                        self._add_line(line.ljust(80) + "│")

                self._add_line("│  └─".ljust(80) + "│")

            # Conditionals
            if self.parser.conditionals:
                self._add_line("│".ljust(80) + "│")
                self._add_line("│  Conditionals:".ljust(80) + "│")
                for condition in self.parser.conditionals:
                    line = f"│  ├──[{condition}]──┐"
                    self._add_line(line.ljust(80) + "│")

            self._add_line("├" + "─" * 78 + "┤")

        # Outputs section
        if self.parser.workflow_outputs:
            self._add_line("│ Outputs:".ljust(80) + "│")
            for var_name, type_info in self.parser.workflow_outputs.items():
                type_str = type_info[0] if isinstance(type_info, tuple) else type_info
                line = f"│   • {var_name}: {type_str}"
                self._add_line(line.ljust(80) + "│")
            self._add_line("├" + "─" * 78 + "┤")

        # Verbose details
        if self.verbose and self.parser.tasks:
            self._add_line("│ Task Details:".ljust(80) + "│")
            for task_name, task_info in self.parser.tasks.items():
                self._add_line(f"│   Task: {task_name}".ljust(80) + "│")
                if task_info.get('inputs'):
                    self._add_line("│     Inputs:".ljust(80) + "│")
                    for var_name, type_info in task_info['inputs'].items():
                        if isinstance(type_info, tuple):
                            type_str, default_val = type_info
                            opt = "?" if "?" in type_str else ""
                            # Mark parameters with defaults as optional
                            if default_val and not opt:
                                opt = "?"
                            if default_val:
                                line = f"│       {var_name}: {type_str}{opt} = {default_val}"
                            else:
                                line = f"│       {var_name}: {type_str}{opt}"
                        else:
                            line = f"│       {var_name}: {type_info}"
                        self._add_line(line.ljust(80) + "│")
                if task_info.get('outputs'):
                    self._add_line("│     Outputs:".ljust(80) + "│")
                    for var_name, type_info in task_info['outputs'].items():
                        type_str = type_info[0] if isinstance(type_info, tuple) else type_info
                        line = f"│       {var_name}: {type_str}"
                        self._add_line(line.ljust(80) + "│")
            self._add_line("├" + "─" * 78 + "┤")

        # Footer
        self._add_line("└" + "─" * 78 + "┘")

        return "\n".join(self.lines)

    def _add_line(self, line: str):
        """Add a line to the output, ensuring it fits width."""
        max_width = 80
        if len(line) > max_width:
            line = line[:max_width-1] + "│"
        self.lines.append(line)


def inspect_workflow(workflow_path: str, verbose: bool = False) -> str:
    """
    Inspect a workflow file and return ASCII diagram.

    Args:
        workflow_path: Path to .wdl file or strategy name (fast/balanced/sensitive)
        verbose: Show detailed parameter information

    Returns:
        ASCII diagram as string
    """
    # Map strategy names to workflow paths
    strategy_map = {
        'fast': 'torchbase/workflows/builtin/fast_typing.wdl',
        'balanced': 'torchbase/workflows/builtin/balanced_typing.wdl',
        'sensitive': 'torchbase/workflows/builtin/sensitive_typing.wdl',
    }

    # Resolve path
    if workflow_path in strategy_map:
        wdl_path = Path(strategy_map[workflow_path])
    else:
        wdl_path = Path(workflow_path)

        # Check if it's a torch directory with main.wdl
        if wdl_path.is_dir():
            main_wdl = wdl_path / 'main.wdl'
            if main_wdl.exists():
                wdl_path = main_wdl
            else:
                raise FileNotFoundError(
                    f"Torch directory has no main.wdl: {wdl_path}"
                )

    # Verify file exists
    if not wdl_path.exists():
        raise FileNotFoundError(f"Workflow not found: {wdl_path}")

    # Read WDL content
    try:
        with open(wdl_path, 'r') as f:
            wdl_content = f.read()
    except Exception as e:
        raise IOError(f"Failed to read workflow: {e}")

    # Parse WDL
    try:
        parser = WDLParser(wdl_content)
    except ValueError as e:
        raise ValueError(f"Failed to parse WDL: {e}")

    # Render diagram
    renderer = WorkflowDiagramRenderer(parser, verbose=verbose)
    return renderer.render()
