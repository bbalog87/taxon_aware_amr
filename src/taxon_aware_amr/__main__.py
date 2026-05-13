"""Allow ``python -m taxon_aware_amr ...`` invocation."""
from .cli import main
raise SystemExit(main())
