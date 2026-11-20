"""Every helper script referenced by the modules exists in scripts/."""
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


class ScriptPathTests(unittest.TestCase):
    def test_module_script_references_exist(self):
        missing = []
        for module in sorted((ROOT / 'modules').glob('*.sh')):
            text = module.read_text()
            self.assertNotRegex(text, r'\$(CONFIG_DIR|SCRIPT_DIR)/scripts/', module.name)
            for name in re.findall(r'/\.\./scripts/([A-Za-z0-9_]+\.py)', text):
                if not (ROOT / 'scripts' / name).is_file():
                    missing.append(f'{module.name}: {name}')
        self.assertEqual(missing, [])


if __name__ == '__main__':
    unittest.main()
