FROM python:3.11-slim

# Instalar gcc, make e outras ferramentas de compilação necessárias para Multigrid
RUN apt-get update && apt-get install -y gcc make && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

# Render injeta a variável PORT, use gunicorn para iniciar
CMD gunicorn --bind 0.0.0.0:$PORT app:app