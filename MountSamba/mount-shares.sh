#!/usr/bin/env bash
MNTDIR="$HOME/mnt"
mkdir -p "$MNTDIR"

CREDFILE="$HOME/.secrets/storage_server_creds"

if [ ! -r "$CREDFILE" ]; then
  echo "ERROR: cannot read creds file at $CREDFILE" >&2
  exit 1
fi

SERVER=$(grep -E '^server'   "$CREDFILE" | cut -d'=' -f2 | xargs)
DOMAIN=$(grep -E '^domain'   "$CREDFILE" | cut -d'=' -f2 | xargs)
USER=$(grep -E '^username'   "$CREDFILE" | cut -d'=' -f2 | xargs)
PASS=$(grep -E '^password'   "$CREDFILE" | cut -d'=' -f2 | xargs)
SHARES=$(grep -E '^shares'   "$CREDFILE" | cut -d'=' -f2 | xargs)

for SHARE in $SHARES; do
  echo "Mounting $SHARE..."
  gio mount "smb://$DOMAIN;$USER:$PASS@$SERVER/$SHARE"
done

echo -e "\nCreating symbolic links in $MNTDIR:"
echo "Might say failed to access share, if the symbolic link already exists .."

for SHARE in $SHARES; do
  DEST="$MNTDIR/$SHARE"

  # Dynamically find gvfs mount regardless of path format
  SRC=$(ls "$XDG_RUNTIME_DIR/gvfs/" 2>/dev/null \
        | grep -i "share=${SHARE,,}" \
        | head -1)

  if [ -z "$SRC" ]; then
    SRC=$(ls "$XDG_RUNTIME_DIR/gvfs/" 2>/dev/null \
          | grep -i "share=$SHARE" \
          | head -1)
  fi

  if [ -n "$SRC" ]; then
    FULL_SRC="$XDG_RUNTIME_DIR/gvfs/$SRC"
    ln -sfn "$FULL_SRC" "$DEST"
    echo "  Linked: $DEST -> $FULL_SRC"
  else
    echo "  ERROR: Could not find gvfs mount for $SHARE"
    echo "  Available gvfs mounts:"
    ls "$XDG_RUNTIME_DIR/gvfs/" 2>/dev/null || echo "    (none)"
  fi
done

echo ""
echo "Active mounts:"
ls -la "$MNTDIR"

echo ""
echo "Verifying mount contents:"
for SHARE in $SHARES; do
  echo "  $MNTDIR/$SHARE:"
  ls "$MNTDIR/$SHARE" 2>/dev/null | head -5 || echo "    (empty or not mounted)"
done
